import os
import pytest
import numpy as np
from qupled.qupled import Qstls as QstlsNative
from qupled.util import Hdf
from qupled.quantum import Qstls


@pytest.fixture
def qstls():
    return Qstls()


@pytest.fixture
def qstls_input():
    return Qstls.Input(1.0, 1.0)


def test_default(qstls):
    assert qstls.hdfFileName is None


def test_compute(qstls, qstls_input, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.Qstls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.quantum.Qstls._checkStatusAndClean")
    mockSave = mocker.patch("qupled.quantum.Qstls._save")
    qstls.compute(qstls_input)
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1


def test_save(qstls, qstls_input, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        scheme = QstlsNative(qstls_input.toNative())
        qstls.hdfFileName = qstls._getHdfFile(scheme.inputs)
        qstls._save(scheme)
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(qstls.hdfFileName)
        inspectData = Hdf().inspect(qstls.hdfFileName)
        expectedEntries = [
            "coupling",
            "degeneracy",
            "theory",
            "error",
            "resolution",
            "cutoff",
            "matsubara",
            "adr",
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
        os.remove(qstls.hdfFileName)


def test_getInitialGuess(mocker):
    arr = np.ones(10)
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read",
        return_value={"wvg": arr, "ssf": arr, "adr": arr, "matsubara": 10},
    )
    guess = Qstls.getInitialGuess("dummyFileName")
    assert np.array_equal(guess.wvg, arr)
    assert np.array_equal(guess.ssf, arr)
    assert np.array_equal(guess.adr, arr)
    assert np.array_equal(guess.matsubara, 10)
