import os
import pytest
import numpy as np
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import Stls
from qupled.classic import ClassicScheme


@pytest.fixture
def stls_instance():
    return Stls(Stls.Input(1.0, 1.0))


def test_default(stls_instance):
    assert issubclass(Stls, ClassicScheme)
    assert issubclass(Stls, qp.Stls)
    assert stls_instance.hdfFileName == "rs1.000_theta1.000_STLS.h5"


def test_compute(stls_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCompute = mocker.patch("qupled.qupled.Stls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.Stls._checkStatusAndClean")
    mockSave = mocker.patch("qupled.classic.Stls._save")
    stls_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSave.call_count == 1


def test_save(stls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    try:
        stls_instance._save()
        assert mockMPIIsRoot.call_count == 2
        assert os.path.isfile(stls_instance.hdfFileName)
        inspectData = Hdf().inspect(stls_instance.hdfFileName)
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
        os.remove(stls_instance.hdfFileName)


def test_getInitialGuess(stls_instance, mocker):
    arr = np.ones(10)
    mockHdfRead = mocker.patch(
        "qupled.util.Hdf.read", return_value={"wvg": arr, "slfc": arr}
    )
    guess = stls_instance.getInitialGuess("dummyFileName")
    assert np.array_equal(guess.wvg, arr)
    assert np.array_equal(guess.slfc, arr)
