import os
import pytest
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import StlsIet, IterativeScheme

@pytest.fixture
def stls_iet_instance():
    return StlsIet(StlsIet.Input(1.0, 1.0, "STLS-HNC"))


def test_default(stls_iet_instance):
    assert issubclass(StlsIet, IterativeScheme)
    assert issubclass(StlsIet, qp.Stls)
    assert stls_iet_instance.hdfFileName == "rs1.000_theta1.000_STLS-HNC.h5"


def test_compute(stls_iet_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.Stls.compute")
    mockCheckStatusAndClean = mocker.patch(
        "qupled.classic.StlsIet._checkStatusAndClean"
    )
    mockSave = mocker.patch("qupled.classic.StlsIet._save")
    stls_iet_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1


def test_save(stls_iet_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        stls_iet_instance._save()
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(stls_iet_instance.hdfFileName)
        inspectData = Hdf().inspect(stls_iet_instance.hdfFileName)
        expectedEntries = [
            "coupling",
            "degeneracy",
            "theory",
            "error",
            "resolution",
            "cutoff",
            "matsubara",
            "bf",
            "idr",
            "sdr",
            "slfc",
            "ssf",
            "ssfHF",
            "wvg",
        ]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(stls_iet_instance.hdfFileName)
