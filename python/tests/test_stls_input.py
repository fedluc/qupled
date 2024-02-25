import os
import pytest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import set_path
import qupled.qupled as qp

@pytest.fixture
def stls_input_instance():
    return qp.StlsInput()


def test_init(stls_input_instance):
    assert issubclass(qp.StlsInput, qp.RpaInput)
    assert hasattr(stls_input_instance, "error")
    assert hasattr(stls_input_instance, "mixing")
    assert hasattr(stls_input_instance, "iet")
    assert hasattr(stls_input_instance, "iterations")
    assert hasattr(stls_input_instance, "outputFrequency")
    assert hasattr(stls_input_instance, "recoveryFile")
    assert hasattr(stls_input_instance, "guess")
    assert hasattr(stls_input_instance.guess, "wvg")
    assert hasattr(stls_input_instance.guess, "slfc")
        
def test_defaults(stls_input_instance):
    assert stls_input_instance.error == 0
    assert stls_input_instance.mixing == 0
    assert stls_input_instance.iet == ""
    assert stls_input_instance.iterations == 0
    assert stls_input_instance.outputFrequency == 0
    assert stls_input_instance.recoveryFile == ""
    assert stls_input_instance.guess.wvg.size == 0
    assert stls_input_instance.guess.slfc.size == 0
    
def test_error(stls_input_instance):
    stls_input_instance.error = 0.001
    error = stls_input_instance.error
    assert error == 0.001
    with pytest.raises(RuntimeError) as excinfo:
        stls_input_instance.error = -0.1
    assert excinfo.value.args[0] == "The minimum error for convergence must be larger than zero"    

def test_mixing(stls_input_instance):
    stls_input_instance.mixing = 0.5
    mixing = stls_input_instance.mixing
    assert mixing == 0.5
    for mixing in [-1, 2]:
        with pytest.raises(RuntimeError) as excinfo:
            stls_input_instance.mixing = -1.0
        assert excinfo.value.args[0] == "The mixing parameter must be a number between zero and one"    

def test_iet(stls_input_instance):
    allowedIet = ["standard", "sqrt", "linear"]
    for iet in allowedIet:
        stls_input_instance.iet = iet
        thisIet = stls_input_instance.iet
        assert thisIet == iet
    with pytest.raises(RuntimeError) as excinfo:
        stls_input_instance.iet = "dummyIET"
    assert excinfo.value.args[0] == "Unknown IET mapping: dummyIET"    

def test_iterations(stls_input_instance):
    stls_input_instance.iterations = 1
    iterations = stls_input_instance.iterations
    assert iterations == 1
    with pytest.raises(RuntimeError) as excinfo:
        stls_input_instance.iterations = -2
    assert excinfo.value.args[0] == "The maximum number of iterations can't be negative"    

def test_outputFrequency(stls_input_instance):
    stls_input_instance.outputFrequency = 1
    outputFrequency = stls_input_instance.outputFrequency
    assert outputFrequency == 1
    with pytest.raises(RuntimeError) as excinfo:
        stls_input_instance.outputFrequency = -3
    assert excinfo.value.args[0] == "The output frequency can't be negative"    

def test_recoveryFile(stls_input_instance):
    stls_input_instance.recoveryFile = "dummyFile"
    recoveryFile = stls_input_instance.recoveryFile
    assert recoveryFile == "dummyFile"

def test_guess(stls_input_instance):
    arr = np.zeros(10)
    guess = qp.SlfcGuess()
    guess.wvg = arr
    guess.slfc = arr
    stls_input_instance.guess = guess
    assert np.array_equal(arr, stls_input_instance.guess.wvg)
    assert np.array_equal(arr, stls_input_instance.guess.slfc)
    with pytest.raises(RuntimeError) as excinfo:
        arr = np.zeros(2)
        guess = qp.SlfcGuess()
        guess.wvg = arr
        guess.slfc = arr
        stls_input_instance.guess = guess
    assert excinfo.value.args[0] == "The initial guess does not contain enough points"
    with pytest.raises(RuntimeError) as excinfo:
        arr1 = np.zeros(10)
        arr2 = np.zeros(11)
        guess = qp.SlfcGuess()
        guess.wvg = arr1
        guess.slfc = arr2
        stls_input_instance.guess = guess
    assert excinfo.value.args[0] == "The initial guess is inconsistent"

def test_isEqual(stls_input_instance):
    thisStls = qp.StlsInput()
    assert stls_input_instance.isEqual(thisStls)
    thisStls.coupling = 2.0
    thisStls.theory = "QSTLS"
    assert not stls_input_instance.isEqual(thisStls)
    
def test_print(stls_input_instance, capfd):
    stls_input_instance.print()
    captured = capfd.readouterr().out
    captured = captured.split("\n")
    assert "Coupling parameter = 0" in captured
    assert "Degeneracy parameter = 0" in captured
    assert "Number of OMP threads = 0" in captured
    assert "Scheme for 2D integrals = " in captured
    assert "Integral relative error = 0" in captured
    assert "Theory to be solved = " in captured
    assert "Guess for chemical potential = 0,0" in captured
    assert "Number of Matsubara frequencies = 0" in captured
    assert "Wave-vector resolution = 0" in captured
    assert "Wave-vector cutoff = 0" in captured
    assert "Iet mapping scheme = " in captured;
    assert "Maximum number of iterations = 0" in captured;
    assert "Minimum error for convergence = 0" in captured;
    assert "Mixing parameter = 0" in captured;
    assert "Output frequency = 0" in captured;
    assert "File with recovery data = " in captured;
