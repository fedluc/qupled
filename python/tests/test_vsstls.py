import os
import pytest
import numpy as np
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import VSStls


@pytest.fixture
def vsstls_instance():
    return VSStls(1.0, 1.0)


def test_default(vsstls_instance):
    assert all(x == y for x, y in zip(vsstls_instance.allowedTheories, ["VSSTLS"]))
    assert vsstls_instance.inputs.coupling == 1.0
    assert vsstls_instance.inputs.degeneracy == 1.0
    assert vsstls_instance.inputs.theory == "VSSTLS"
    assert all(x == y for x, y in zip(vsstls_instance.inputs.chemicalPotential, [-100.0, 100.0]))
    assert vsstls_instance.inputs.cutoff == 10.0
    assert vsstls_instance.inputs.error == 1.0e-5
    assert vsstls_instance.inputs.mixing == 1.0
    assert vsstls_instance.inputs.iterations == 1000
    assert vsstls_instance.inputs.matsubara == 128
    assert vsstls_instance.inputs.outputFrequency == 10
    assert vsstls_instance.inputs.recoveryFile == ""
    assert vsstls_instance.inputs.resolution == 0.1
    assert all(x == y for x, y in zip(vsstls_instance.inputs.alpha, [0.5, 1.0]))
    assert vsstls_instance.inputs.couplingResolution == 0.01
    assert vsstls_instance.inputs.degeneracyResolution == 0.01
    assert vsstls_instance.inputs.errorAlpha == 1.0e-3
    assert vsstls_instance.inputs.iterationsAlpha == 50
    assert vsstls_instance.inputs.intError == 1.0e-5
    assert vsstls_instance.inputs.threads == 1
    assert vsstls_instance.scheme is None
    assert vsstls_instance.hdfFileName is None

                                      
def test_set_input():
    vsstls = VSStls(2.0, 0.5, [-5, 5], 20, 1.0e-8, 0.5, 100,
                    32, 100, "recoveryFile", 0.01, [0.8, 1.2],
                    0.1, 0.2, 1.0e-5, 100, 1.0e-8, 9)
    assert vsstls.inputs.coupling == 2.0
    assert vsstls.inputs.degeneracy == 0.5
    assert vsstls.inputs.theory == "VSSTLS"
    assert all(x == y for x, y in zip(vsstls.inputs.chemicalPotential, [-5.0, 5.0]))
    assert vsstls.inputs.cutoff == 20
    assert vsstls.inputs.error == 1.0e-8
    assert vsstls.inputs.mixing == 0.5
    assert vsstls.inputs.iterations == 100
    assert vsstls.inputs.matsubara == 32
    assert vsstls.inputs.outputFrequency == 100
    assert vsstls.inputs.recoveryFile == "recoveryFile"
    assert vsstls.inputs.resolution == 0.01
    assert all(x == y for x, y in zip(vsstls.inputs.alpha, [0.8, 1.2]))
    assert vsstls.inputs.couplingResolution == 0.1
    assert vsstls.inputs.degeneracyResolution == 0.2
    assert vsstls.inputs.errorAlpha == 1.0e-5
    assert vsstls.inputs.iterationsAlpha == 100
    assert vsstls.inputs.intError == 1.0e-8
    assert vsstls.inputs.threads == 9
    assert vsstls.scheme is None
    assert vsstls.hdfFileName is None
    
def test_compute(vsstls_instance, mocker):
    mockMPITime = mocker.patch("qupled.util.MPI.timer", return_value=0)
    mockMPIBarrier = mocker.patch("qupled.util.MPI.barrier")
    mockCheckInputs = mocker.patch("qupled.classic.VSStls._checkInputs")
    mockCompute = mocker.patch("qupled.qupled.VSStls.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.VSStls._checkStatusAndClean")
    mockSetHdfFile = mocker.patch("qupled.classic.VSStls._setHdfFile")
    mockSave = mocker.patch("qupled.classic.VSStls._save")
    vsstls_instance.compute()
    assert mockMPITime.call_count == 2
    assert mockMPIBarrier.call_count == 1
    assert mockCheckInputs.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1

    
def test_save(vsstls_instance, mocker):
    mockMPIIsRoot = mocker.patch("qupled.util.MPI.isRoot")
    vsstls_instance.scheme = qp.VSStls(vsstls_instance.inputs)
    vsstls_instance._setHdfFile()
    try:
        vsstls_instance._save()
        assert mockMPIIsRoot.call_count == 3
        assert os.path.isfile(vsstls_instance.hdfFileName)
        inspectData = Hdf().inspect(vsstls_instance.hdfFileName)
        expectedEntries = ["coupling", "degeneracy", "theory", "error",
                           "resolution", "cutoff", "matsubara",
                           "idr", "sdr", "slfc",  "ssf", "ssfHF",
                           "wvg", "fxcGrid", "fxci"]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(vsstls_instance.hdfFileName)


def test_setFreeEnergyIntegrand(vsstls_instance, mocker):
    arr1D = np.ones(10)
    arr2D = np.ones((3, 10))
    mockHdfRead = mocker.patch("qupled.util.Hdf.read", return_value={"fxcGrid" : arr1D, "fxci" : arr2D})
    vsstls_instance.setFreeEnergyIntegrand("dummyFileName")
    assert np.array_equal(vsstls_instance.inputs.freeEnergyIntegrand.grid, arr1D)
    assert np.array_equal(vsstls_instance.inputs.freeEnergyIntegrand.integrand, arr2D)
