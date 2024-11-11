import os
import pytest
import numpy as np
import qupled.qupled as qp


@pytest.fixture
def qstls_input_instance():
    return qp.QstlsInput()


def test_init(qstls_input_instance):
    assert issubclass(qp.QstlsInput, qp.RpaInput)
    assert hasattr(qstls_input_instance, "guess")
    assert hasattr(qstls_input_instance.guess, "wvg")
    assert hasattr(qstls_input_instance.guess, "ssf")
    assert hasattr(qstls_input_instance.guess, "adr")
    assert hasattr(qstls_input_instance.guess, "matsubara")
    assert hasattr(qstls_input_instance, "fixed")
    assert hasattr(qstls_input_instance, "fixediet")


def test_defaults(qstls_input_instance):
    assert qstls_input_instance.guess.wvg.size == 0
    assert qstls_input_instance.guess.ssf.size == 0
    assert qstls_input_instance.guess.adr.size == 0
    assert qstls_input_instance.guess.matsubara == 0
    assert qstls_input_instance.fixed == ""
    assert qstls_input_instance.fixediet == ""


def test_fixed(qstls_input_instance):
    qstls_input_instance.fixed = "fixedFile"
    fixed = qstls_input_instance.fixed
    assert fixed == "fixedFile"


def test_fixed(qstls_input_instance):
    qstls_input_instance.fixediet = "fixedFile"
    fixed = qstls_input_instance.fixediet
    assert fixed == "fixedFile"


def test_guess(qstls_input_instance):
    arr = np.zeros(10)
    guess = qp.QstlsGuess()
    guess.wvg = arr
    guess.ssf = arr
    qstls_input_instance.guess = guess
    assert np.array_equal(arr, qstls_input_instance.guess.wvg)
    assert np.array_equal(arr, qstls_input_instance.guess.ssf)
    with pytest.raises(RuntimeError) as excinfo:
        arr1 = np.zeros(10)
        arr2 = np.zeros(11)
        guess = qp.QstlsGuess()
        guess.wvg = arr1
        guess.ssf = arr2
        qstls_input_instance.guess = guess
    assert excinfo.value.args[0] == "The initial guess is inconsistent"


def test_guessIet(qstls_input_instance):
    arr1 = np.zeros(10)
    arr2 = np.zeros((10, 4))
    matsubara = 4
    guess = qp.QstlsGuess()
    guess.wvg = arr1
    guess.ssf = arr1
    guess.adr = arr2
    guess.matsubara = matsubara
    qstls_input_instance.guess = guess
    assert np.array_equal(arr1, qstls_input_instance.guess.wvg)
    assert np.array_equal(arr1, qstls_input_instance.guess.ssf)
    assert np.array_equal(arr2, qstls_input_instance.guess.adr)
    assert guess.matsubara == matsubara
    for arr2 in [np.zeros((11, 4)), np.zeros((10, 5))]:
        with pytest.raises(RuntimeError) as excinfo:
            guess.adr = arr2
            qstls_input_instance.guess = guess
        assert excinfo.value.args[0] == "The initial guess is inconsistent"


def test_isEqual_default(qstls_input_instance):
    thisQstls = qp.QstlsInput()
    assert not qstls_input_instance.isEqual(thisQstls)


def test_isEqual_nonDefault(qstls_input_instance):
    thisQstls = qp.QstlsInput()
    thisQstls.coupling = 2.0
    thisQstls.degeneracy = 1.0
    thisQstls.intError = 0.1
    thisQstls.threads = 1
    thisQstls.theory = "STLS"
    thisQstls.matsubara = 1
    thisQstls.resolution = 0.1
    thisQstls.cutoff = 1.0
    thisQstls.error = 0.1
    thisQstls.mixing = 1.0
    thisQstls.outputFrequency = 1
    assert thisQstls.isEqual(thisQstls)


def test_print(qstls_input_instance, capfd):
    qstls_input_instance.print()
    captured = capfd.readouterr().out
    captured = captured.split("\n")
    assert "Coupling parameter = nan" in captured
    assert "Degeneracy parameter = nan" in captured
    assert "Number of OMP threads = 0" in captured
    assert "Scheme for 2D integrals = " in captured
    assert "Integral relative error = nan" in captured
    assert "Theory to be solved = " in captured
    assert "Guess for chemical potential = " in captured
    assert "Number of Matsubara frequencies = 0" in captured
    assert "Wave-vector resolution = nan" in captured
    assert "Wave-vector cutoff = nan" in captured
    assert "Iet mapping scheme = " in captured
    assert "Maximum number of iterations = 0" in captured
    assert "Minimum error for convergence = nan" in captured
    assert "Mixing parameter = nan" in captured
    assert "Output frequency = 0" in captured
    assert "File with recovery data = " in captured
    assert "File with fixed adr component = " in captured
    assert "File with fixed adr component (iet) = " in captured
