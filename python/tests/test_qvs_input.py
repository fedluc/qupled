import os
import pytest
import numpy as np
import qupled.qupled as qp


@pytest.fixture
def qvsstls_input_instance():
    return qp.QVSStlsInput()


def test_init(qvsstls_input_instance):
    assert issubclass(qp.QVSStlsInput, qp.VSInput)
    assert hasattr(qvsstls_input_instance, "errorAlpha")
    assert hasattr(qvsstls_input_instance, "iterationsAlpha")
    assert hasattr(qvsstls_input_instance, "alpha")
    assert hasattr(qvsstls_input_instance, "couplingResolution")
    assert hasattr(qvsstls_input_instance, "degeneracyResolution")
    assert hasattr(qvsstls_input_instance, "freeEnergyIntegrand")
    assert hasattr(qvsstls_input_instance, "guess")
    assert hasattr(qvsstls_input_instance.guess, "wvg")
    assert hasattr(qvsstls_input_instance.guess, "ssf")
    assert hasattr(qvsstls_input_instance.guess, "adr")
    assert hasattr(qvsstls_input_instance.guess, "matsubara")
    assert hasattr(qvsstls_input_instance, "fixed")
    assert hasattr(qvsstls_input_instance.freeEnergyIntegrand, "grid")
    assert hasattr(qvsstls_input_instance.freeEnergyIntegrand, "integrand")


def test_defaults(qvsstls_input_instance):
    assert np.isnan(qvsstls_input_instance.errorAlpha)
    assert qvsstls_input_instance.iterationsAlpha == 0
    assert qvsstls_input_instance.alpha.size == 0
    assert np.isnan(qvsstls_input_instance.couplingResolution)
    assert np.isnan(qvsstls_input_instance.degeneracyResolution)
    assert qvsstls_input_instance.freeEnergyIntegrand.grid.size == 0
    assert qvsstls_input_instance.freeEnergyIntegrand.integrand.size == 0
    assert qvsstls_input_instance.guess.wvg.size == 0
    assert qvsstls_input_instance.guess.ssf.size == 0
    assert qvsstls_input_instance.guess.adr.size == 0
    assert qvsstls_input_instance.guess.matsubara == 0
    assert qvsstls_input_instance.fixed == ""


def test_fixed(qvsstls_input_instance):
    qvsstls_input_instance.fixed = "fixedFile"
    fixed = qvsstls_input_instance.fixed
    assert fixed == "fixedFile"


def test_errorAlpha(qvsstls_input_instance):
    qvsstls_input_instance.errorAlpha = 0.001
    errorAlpha = qvsstls_input_instance.errorAlpha
    assert errorAlpha == 0.001
    with pytest.raises(RuntimeError) as excinfo:
        qvsstls_input_instance.errorAlpha = -0.1
    assert (
        excinfo.value.args[0]
        == "The minimum error for convergence must be larger than zero"
    )


def test_iterationsAlpha(qvsstls_input_instance):
    qvsstls_input_instance.iterationsAlpha = 1
    iterationsAlpha = qvsstls_input_instance.iterationsAlpha
    assert iterationsAlpha == 1
    with pytest.raises(RuntimeError) as excinfo:
        qvsstls_input_instance.iterationsAlpha = -2
    assert excinfo.value.args[0] == "The maximum number of iterations can't be negative"


def test_alpha(qvsstls_input_instance):
    qvsstls_input_instance.alpha = [-10, 10]
    alpha = qvsstls_input_instance.alpha
    assert all(x == y for x, y in zip(alpha, [-10, 10]))
    for a in [[-1.0], [1, 2, 3], [10, -10]]:
        with pytest.raises(RuntimeError) as excinfo:
            qvsstls_input_instance.alpha = a
        assert excinfo.value.args[0] == "Invalid guess for free parameter calculation"


def test_couplingResolution(qvsstls_input_instance):
    qvsstls_input_instance.couplingResolution = 0.01
    couplingResolution = qvsstls_input_instance.couplingResolution
    assert couplingResolution == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        qvsstls_input_instance.couplingResolution = -0.1
    assert (
        excinfo.value.args[0]
        == "The coupling parameter resolution must be larger than zero"
    )


def test_degeneracyResolution(qvsstls_input_instance):
    qvsstls_input_instance.degeneracyResolution = 0.01
    degeneracyResolution = qvsstls_input_instance.degeneracyResolution
    assert degeneracyResolution == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        qvsstls_input_instance.degeneracyResolution = -0.1
    assert (
        excinfo.value.args[0]
        == "The degeneracy parameter resolution must be larger than zero"
    )


def test_freeEnergyIntegrand(qvsstls_input_instance):
    arr1 = np.zeros(10)
    arr2 = np.zeros((3, 10))
    fxc = qp.FreeEnergyIntegrand()
    fxc.grid = arr1
    fxc.integrand = arr2
    qvsstls_input_instance.freeEnergyIntegrand = fxc
    assert np.array_equal(arr1, qvsstls_input_instance.freeEnergyIntegrand.grid)
    assert np.array_equal(arr2, qvsstls_input_instance.freeEnergyIntegrand.integrand)


def test_freeEnergyIntegrand_Inconsistent(qvsstls_input_instance):
    with pytest.raises(RuntimeError) as excinfo:
        arr1 = np.zeros(10)
        arr2 = np.zeros((3, 11))
        fxc = qp.FreeEnergyIntegrand()
        fxc.grid = arr1
        fxc.integrand = arr2
        qvsstls_input_instance.freeEnergyIntegrand = fxc
    assert excinfo.value.args[0] == "The free energy integrand is inconsistent"


def test_isEqual_default(qvsstls_input_instance):
    assert not qvsstls_input_instance.isEqual(qvsstls_input_instance)


def test_isEqual(qvsstls_input_instance):
    thisQVSStls = qp.QVSStlsInput()
    thisQVSStls.coupling = 2.0
    thisQVSStls.degeneracy = 1.0
    thisQVSStls.intError = 0.1
    thisQVSStls.threads = 1
    thisQVSStls.theory = "STLS"
    thisQVSStls.matsubara = 1
    thisQVSStls.resolution = 0.1
    thisQVSStls.cutoff = 1.0
    thisQVSStls.error = 0.1
    thisQVSStls.mixing = 1.0
    thisQVSStls.outputFrequency = 1
    thisQVSStls.couplingResolution = 0.1
    thisQVSStls.degeneracyResolution = 0.1
    thisQVSStls.errorAlpha = 0.1
    thisQVSStls.iterationsAlpha = 1
    assert thisQVSStls.isEqual(thisQVSStls)


def test_print(qvsstls_input_instance, capfd):
    qvsstls_input_instance.print()
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
    assert "Guess for the free parameter = " in captured
    assert "Resolution for the coupling parameter grid = nan" in captured
    assert "Resolution for the degeneracy parameter grid = nan" in captured
    assert "Minimum error for convergence (alpha) = nan" in captured
    assert "Maximum number of iterations (alpha) = 0" in captured
    assert "File with fixed adr component = " in captured
