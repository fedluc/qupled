import os
import pytest
import numpy as np
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import Stls


@pytest.fixture
def stls_instance():
    return Stls(1.0, 1.0)


def test_default(stls_instance):
    assert all(x == y for x, y in zip(stls_instance.allowedTheories, ["STLS"]))
    assert stls_instance.inputs.coupling == 1.0
    assert stls_instance.inputs.degeneracy == 1.0
    assert stls_instance.inputs.theory == "STLS"
    assert all(x == y for x, y in zip(stls_instance.inputs.chemicalPotential, [-10.0, 10.0]))
    assert stls_instance.inputs.cutoff == 10.0
    assert stls_instance.inputs.error == 1.0e-5
    assert stls_instance.inputs.mixing == 1.0
    assert stls_instance.inputs.iterations == 1000
    assert stls_instance.inputs.matsubara == 128
    assert stls_instance.inputs.outputFrequency == 10
    assert stls_instance.inputs.recoveryFile == ""
    assert stls_instance.inputs.resolution == 0.1
    assert stls_instance.inputs.intError == 1.0e-5
    assert stls_instance.inputs.threads == 1
    assert stls_instance.scheme is None
    assert stls_instance.hdfFileName is None

                                      
def test_set_input():
    stls = Stls(2.0, 0.5, [-5, 5], 20, 1.0e-8, 0.5, None,
                100, 32, 100, "recoveryFile", 0.01)
    assert stls.inputs.coupling == 2.0
    assert stls.inputs.degeneracy == 0.5
    assert stls.inputs.theory == "STLS"
    assert all(x == y for x, y in zip(stls.inputs.chemicalPotential, [-5.0, 5.0]))
    assert stls.inputs.cutoff == 20
    assert stls.inputs.error == 1.0e-8
    assert stls.inputs.mixing == 0.5
    assert stls.inputs.iterations == 100
    assert stls.inputs.matsubara == 32
    assert stls.inputs.outputFrequency == 100
    assert stls.inputs.recoveryFile == "recoveryFile"
    assert stls.inputs.resolution == 0.01
    
    
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
    stls_instance.scheme = qp.Stls(stls_instance.inputs)
    stls_instance._setHdfFile()
    try:
        stls_instance._save()
        assert mockMPIIsRoot.call_count == 2
        assert os.path.isfile(stls_instance.hdfFileName)
        inspectData = Hdf().inspect(stls_instance.hdfFileName)
        expectedEntries = ["coupling", "degeneracy", "theory", "error",
                           "resolution", "cutoff", "matsubara",
                           "idr", "sdr", "slfc",
                           "ssf", "ssfHF", "wvg"]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(stls_instance.hdfFileName)


def test_setGuess(stls_instance, mocker):
    arr = np.ones(10)
    mockHdfRead = mocker.patch("qupled.util.Hdf.read", return_value={"wvg" : arr, "slfc" : arr})
    stls_instance.setGuess("dummyFileName")
    assert np.array_equal(stls_instance.inputs.guess.wvg, arr)
    assert np.array_equal(stls_instance.inputs.guess.slfc, arr)
