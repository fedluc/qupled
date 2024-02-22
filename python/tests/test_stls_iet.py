import os
import pytest
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import StlsIet


@pytest.fixture
def stls_iet_instance():
    return StlsIet(1.0, 1.0, "STLS-HNC")


def test_default(stls_iet_instance):
    assert all(x == y for x, y in zip(stls_iet_instance.allowedTheories, ["STLS-HNC",
                                                                          "STLS-IOI",
                                                                          "STLS-LCT"]))
    assert stls_iet_instance.inputs.coupling == 1.0
    assert stls_iet_instance.inputs.degeneracy == 1.0
    assert stls_iet_instance.inputs.theory == "STLS-HNC"
    assert all(x == y for x, y in zip(stls_iet_instance.inputs.chemicalPotential, [-10.0, 10.0]))
    assert stls_iet_instance.inputs.cutoff == 10.0
    assert stls_iet_instance.inputs.error == 1.0e-5
    assert stls_iet_instance.inputs.iet == "standard"
    assert stls_iet_instance.inputs.mixing == 1.0
    assert stls_iet_instance.inputs.iterations == 1000    
    assert stls_iet_instance.inputs.matsubara == 128
    assert stls_iet_instance.inputs.outputFrequency == 10
    assert stls_iet_instance.inputs.recoveryFile == ""
    assert stls_iet_instance.inputs.resolution == 0.1
    assert stls_iet_instance.inputs.int2DScheme == "full"
    assert stls_iet_instance.inputs.intError == 1.0e-5
    assert stls_iet_instance.inputs.threads == 1
    assert stls_iet_instance.scheme is None
    assert stls_iet_instance.hdfFileName is None

                                      
def test_set_input():
    stls_iet = StlsIet(2.0, 0.5, "STLS-LCT", [-5, 5], 20, 1.0e-8, "sqrt",
                       0.5, None, 100, 32, 100, "recoveryFile", 0.01,
                       "segregated")
    assert stls_iet.inputs.coupling == 2.0
    assert stls_iet.inputs.degeneracy == 0.5
    assert stls_iet.inputs.theory == "STLS-LCT"
    assert all(x == y for x, y in zip(stls_iet.inputs.chemicalPotential, [-5.0, 5.0]))
    assert stls_iet.inputs.cutoff == 20.0
    assert stls_iet.inputs.error == 1.0e-8
    assert stls_iet.inputs.iet == "sqrt"
    assert stls_iet.inputs.mixing == 0.5
    assert stls_iet.inputs.iterations == 100    
    assert stls_iet.inputs.matsubara == 32
    assert stls_iet.inputs.outputFrequency == 100
    assert stls_iet.inputs.recoveryFile == "recoveryFile"
    assert stls_iet.inputs.resolution == 0.01
    assert stls_iet.inputs.int2DScheme == "segregated"
    
    
def test_compute(stls_iet_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.classic.StlsIet._checkInputs")
    mockCompute = mocker.patch("qupled.qupled.Stls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.StlsIet._checkStatusAndClean")
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
        expectedEntries = ["coupling", "degeneracy", "theory", "error",
                           "resolution", "cutoff", "matsubara", "bf",
                           "idr", "sdr", "slfc",
                           "ssf", "ssfHF", "wvg"]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(stls_iet_instance.hdfFileName)
