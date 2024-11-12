import os
import pytest
from qupled.qupled import Stls as StlsNative
from qupled.util import Hdf
from qupled.classic import StlsIet


@pytest.fixture
def stls_iet():
    return StlsIet()


@pytest.fixture
def stls_iet_input():
    return StlsIet.Input(1.0, 1.0, "STLS-HNC")


def test_default(stls_iet):
    assert stls_iet.hdfFileName is None


def test_compute(stls_iet, stls_iet_input, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.Stls.compute")
    mockCheckStatusAndClean = mocker.patch(
        "qupled.classic.StlsIet._checkStatusAndClean"
    )
    mockSave = mocker.patch("qupled.classic.StlsIet._save")
    stls_iet.compute(stls_iet_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1


def test_save(stls_iet, stls_iet_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        scheme = StlsNative(stls_iet_input.toNative())
        stls_iet.hdfFileName = stls_iet._getHdfFile(scheme.inputs)
        stls_iet._save(scheme)
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(stls_iet.hdfFileName)
        inspectData = Hdf().inspect(stls_iet.hdfFileName)
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
        os.remove(stls_iet.hdfFileName)
