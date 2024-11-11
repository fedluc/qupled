import os
import pytest
import numpy as np
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.quantum import Qstls, QuantumIterativeScheme


@pytest.fixture
def qstls_instance():
    return Qstls(Qstls.Input(1.0, 1.0))


def test_default(qstls_instance):
    assert issubclass(Qstls, QuantumIterativeScheme)
    assert issubclass(Qstls, qp.Stls)
    assert qstls_instance.hdfFileName == "rs1.000_theta1.000_QSTLS.h5"


def test_compute(qstls_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.Qstls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.quantum.Qstls._checkStatusAndClean")
    mockSave = mocker.patch("qupled.quantum.Qstls._save")
    qstls_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1


def test_save(qstls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        qstls_instance._save()
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(qstls_instance.hdfFileName)
        inspectData = Hdf().inspect(qstls_instance.hdfFileName)
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
        os.remove(qstls_instance.hdfFileName)


def test_getInitialGuess(qstls_instance, mocker):
    arr = np.ones(10)
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read", return_value={"wvg": arr, "ssf": arr}
    )
    guess = Qstls.getInitialGuess("dummyFileName")
    assert np.array_equal(guess.wvg, arr)
    assert np.array_equal(guess.ssf, arr)
