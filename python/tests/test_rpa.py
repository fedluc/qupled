import os
import pytest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import set_path
import qupled.qupled as qp
from qupled.util import Hdf
from qupled.classic import Rpa


@pytest.fixture
def rpa_instance():
    return Rpa(1.0, 1.0)


def test_default(rpa_instance):
    assert all(x == y for x, y in zip(rpa_instance.allowedTheories, ["RPA"]))
    assert rpa_instance.inputs.coupling == 1.0
    assert rpa_instance.inputs.degeneracy == 1.0
    assert rpa_instance.inputs.theory == "RPA"
    assert all(x == y for x, y in zip(rpa_instance.inputs.chemicalPotential, [-10.0, 10.0]))
    assert rpa_instance.inputs.cutoff == 10.0
    assert rpa_instance.inputs.matsubara == 128
    assert rpa_instance.inputs.resolution == 0.1
    assert rpa_instance.inputs.intError == 1.0e-5
    assert rpa_instance.inputs.threads == 1
    assert rpa_instance.scheme is None
    assert rpa_instance.hdfFileName is None

                                      
def test_set_input():
    rpa = Rpa(2.0, 0.5, [-5, 5], 20, 32, 0.01)
    assert rpa.inputs.coupling == 2.0
    assert rpa.inputs.degeneracy == 0.5
    assert rpa.inputs.theory == "RPA"
    assert all(x == y for x, y in zip(rpa.inputs.chemicalPotential, [-5.0, 5.0]))
    assert rpa.inputs.cutoff == 20
    assert rpa.inputs.matsubara == 32
    assert rpa.inputs.resolution == 0.01

    
def test_compute(rpa_instance, mocker):
    mockCheckInputs = mocker.patch("qupled.classic.Rpa._checkInputs")
    mockCompute = mocker.patch("qupled.qupled.Rpa.compute")
    mockCheckStatusAndClean = mocker.patch("qupled.classic.Rpa._checkStatusAndClean")
    mockSetHdfFile = mocker.patch("qupled.classic.Rpa._setHdfFile")
    mockSave = mocker.patch("qupled.classic.Rpa._save")
    rpa_instance.compute()
    assert mockCheckInputs.call_count == 1
    assert mockCompute.call_count == 1
    assert mockCheckStatusAndClean.call_count == 1
    assert mockSetHdfFile.call_count == 1
    assert mockSave.call_count == 1

    
def test_checkStatusAndClean(rpa_instance, mocker, capsys):
    mockCheckInputs = mocker.patch("os.remove")
    rpa_instance.scheme = qp.Rpa(rpa_instance.inputs)
    rpa_instance._checkStatusAndClean(0)
    captured = capsys.readouterr()
    assert "Dielectric theory solved successfully!\n" in captured
    with pytest.raises(SystemExit) as excinfo:
        rpa_instance._checkStatusAndClean(1)
    assert excinfo.value.code == "Error while solving the dielectric theory"


def test_setHdfFile(rpa_instance):
    rpa_instance._setHdfFile()
    assert rpa_instance.hdfFileName == "rs1.000_theta1.000_RPA.h5"

    
def test_save(rpa_instance):
    rpa_instance.scheme = qp.Rpa(rpa_instance.inputs)
    rpa_instance._setHdfFile()
    try:
        rpa_instance._save()
        assert os.path.isfile(rpa_instance.hdfFileName)
        inspectData = Hdf().inspect(rpa_instance.hdfFileName)
        expectedEntries = ["coupling", "degeneracy", "theory", "resolution",
                           "cutoff", "matsubara", "idr", "sdr", "slfc",
                           "ssf", "ssfHF", "wvg"]
        for entry in expectedEntries:
            assert entry in inspectData
    finally:
        os.remove(rpa_instance.hdfFileName)

        
def test_computeRdf(rpa_instance, mocker):
    mockCheckSolution = mocker.patch("qupled.classic.Rpa._checkSolution")
    mockMPIGetRank = mocker.patch("qupled.util.MPI.getRank", return_value=0)
    mockComputeRdf = mocker.patch("qupled.util.Hdf.computeRdf")
    rpa_instance.computeRdf()
    assert mockCheckSolution.call_count == 1
    assert mockMPIGetRank.call_count == 1
    assert mockComputeRdf.call_count == 1


def test_computeInternalEnergy(rpa_instance, mocker):
    mockCheckSolution = mocker.patch("qupled.classic.Rpa._checkSolution")
    mockComputeInternalEnergy = mocker.patch("qupled.qupled.computeInternalEnergy")
    rpa_instance.scheme = qp.Rpa(rpa_instance.inputs)
    rpa_instance.computeInternalEnergy()
    assert mockCheckSolution.call_count == 1
    assert mockComputeInternalEnergy.call_count == 1

def test_plot(rpa_instance, mocker):
    mockCheckSolution = mocker.patch("qupled.classic.Rpa._checkSolution")
    mockComputeRdf = mocker.patch("qupled.classic.Rpa.computeRdf")
    mockPlot = mocker.patch("qupled.util.Hdf.plot")
    rpa_instance.plot(["ssf", "idr"])
    assert mockCheckSolution.call_count == 1
    assert mockComputeRdf.call_count == 0
    assert mockPlot.call_count == 1
    rpa_instance.plot(["ssf", "rdf"])
    assert mockCheckSolution.call_count == 2
    assert mockComputeRdf.call_count == 1
    assert mockPlot.call_count == 2

def test_checkSolution(rpa_instance):
    with pytest.raises(SystemExit) as excinfo:
        rpa_instance._checkSolution("dummy action")
    assert excinfo.value.code == "No solution to dummy action"
    rpa_instance.scheme = qp.Rpa(rpa_instance.inputs)
    try:
        rpa_instance._checkSolution("dummy action")
    except:
        assert False
