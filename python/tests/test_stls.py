import os
import pytest
import numpy as np
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import Stls
from qupled.classic import ClassicSchemeNew


@pytest.fixture
def stls_instance():
    inputs = qp.StlsInput()
    inputs.coupling = 1.0
    inputs.degeneracy = 1.0
    inputs.theory = "STLS"
    inputs.chemicalPotential = [-10, 10]
    inputs.cutoff = 10.0
    inputs.matsubara = 128
    inputs.resolution = 0.1
    inputs.intError = 1.0e-5
    inputs.threads = 1
    inputs.error = 1.0e-5
    inputs.mixing = 1.0
    inputs.iterations = 1000
    inputs.outputFrequency = 10
    return Stls(inputs)


def test_default(stls_instance):
    assert issubclass(Stls, ClassicSchemeNew)
    assert issubclass(Stls, qp.Stls)
    assert all(x == y for x, y in zip(stls_instance.allowedTheories, ["STLS"]))
    assert stls_instance.hdfFileName is None


def test_compute(stls_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.classic.Stls._checkInputs")
    mockCompute = mocker.patch("qupled.qupled.Stls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.Stls._checkStatusAndClean")
    mockSetHdfFile = mocker.patch("qupled.classic.Stls._setHdfFile")
    mockSave = mocker.patch("qupled.classic.Stls._save")
    stls_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCheckInputs.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1


def test_save(stls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    stls_instance._setHdfFile()
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
