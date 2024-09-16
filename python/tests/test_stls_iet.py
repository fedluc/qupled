import os
import pytest
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import StlsIet
from qupled.classic import Stls


@pytest.fixture
def stls_iet_instance():
    return StlsIet(StlsIet.Input(1.0, 1.0, "STLS-HNC"))


def test_default(stls_iet_instance):
    assert issubclass(StlsIet, Stls)
    assert all(
        x == y
        for x, y in zip(
            stls_iet_instance.allowedTheories, ["STLS-HNC", "STLS-IOI", "STLS-LCT"]
        )
    )


def test_compute(stls_iet_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.classic.StlsIet._checkInputs")
    mockCompute = mocker.patch("qupled.qupled.Stls.compute")
    mockCheckStatusAndClean = mocker.patch(
        "qupled.classic.StlsIet._checkStatusAndClean"
    )
    mockSetHdfFile = mocker.patch("qupled.classic.StlsIet._setHdfFile")
    mockSave = mocker.patch("qupled.classic.StlsIet._save")
    stls_iet_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCheckInputs.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1


def test_save(stls_iet_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    stls_iet_instance.scheme = qp.Stls(stls_iet_instance.inputs)
    stls_iet_instance._setHdfFile()
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
