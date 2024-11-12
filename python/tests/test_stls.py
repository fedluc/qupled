import os
import pytest
import numpy as np
from qupled.qupled import Stls as StlsNative
from qupled.util import Hdf
from qupled.classic import Stls


@pytest.fixture
def stls():
    return Stls()


@pytest.fixture
def stls_input():
    return Stls.Input(1.0, 1.0)


def test_default(stls):
    assert stls.hdfFileName is None


def test_compute(stls, stls_input, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.Stls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.Stls._checkStatusAndClean")
    mockSave = mocker.patch("qupled.classic.Stls._save")
    stls.compute(stls_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1


def test_save(stls, stls_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        scheme = StlsNative(stls_input.toNative())
        stls.hdfFileName = stls._getHdfFile(scheme.inputs)
        stls._save(scheme)
        assert mockMPIIsRoot.call_count == 2
        assert os.path.isfile(stls.hdfFileName)
        inspectData = Hdf().inspect(stls.hdfFileName)
        expectedEntries = [
            "coupling",
            "degeneracy",
            "theory",
            "error",
            "resolution",
            "cutoff",
            "matsubara",
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
        os.remove(stls.hdfFileName)


def test_getInitialGuess(mocker):
    arr = np.ones(10)
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read", return_value={"wvg": arr, "slfc": arr}
    )
    guess = Stls.getInitialGuess("dummyFileName")
    assert np.array_equal(guess.wvg, arr)
    assert np.array_equal(guess.slfc, arr)
