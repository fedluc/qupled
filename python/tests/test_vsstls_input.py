import os
import pytest
import numpy as np
import set_path
import qupled.qupled as qp


@pytest.fixture
def vsstls_input_instance():
    return qp.VSStlsInput()


def test_init(vsstls_input_instance):
    assert issubclass(qp.VSStlsInput, qp.VSInput)
    assert hasattr(vsstls_input_instance, "errorAlpha")
    assert hasattr(vsstls_input_instance, "iterationsAlpha")
    assert hasattr(vsstls_input_instance, "alpha")
    assert hasattr(vsstls_input_instance, "couplingResolution")
    assert hasattr(vsstls_input_instance, "degeneracyResolution")
    assert hasattr(vsstls_input_instance, "freeEnergyIntegrand")
    assert hasattr(vsstls_input_instance.freeEnergyIntegrand, "grid")
    assert hasattr(vsstls_input_instance.freeEnergyIntegrand, "integrand")
    assert hasattr(vsstls_input_instance, "iet")
    assert hasattr(vsstls_input_instance, "guess")


def test_defaults(vsstls_input_instance):
    assert vsstls_input_instance.errorAlpha == 0
    assert vsstls_input_instance.iterationsAlpha == 0
    assert all(x == y for x, y in zip(vsstls_input_instance.alpha, [0, 0]))
    assert vsstls_input_instance.couplingResolution == 0
    assert vsstls_input_instance.degeneracyResolution == 0
    assert vsstls_input_instance.freeEnergyIntegrand.grid.size == 0
    assert vsstls_input_instance.freeEnergyIntegrand.integrand.size == 0


def test_errorAlpha(vsstls_input_instance):
    vsstls_input_instance.errorAlpha = 0.001
    errorAlpha = vsstls_input_instance.errorAlpha
    assert errorAlpha == 0.001
    with pytest.raises(RuntimeError) as excinfo:
        vsstls_input_instance.errorAlpha = -0.1
    assert (
        excinfo.value.args[0]
        == "The minimum error for convergence must be larger than zero"
    )


def test_iterationsAlpha(vsstls_input_instance):
    vsstls_input_instance.iterationsAlpha = 1
    iterationsAlpha = vsstls_input_instance.iterationsAlpha
    assert iterationsAlpha == 1
    with pytest.raises(RuntimeError) as excinfo:
        vsstls_input_instance.iterationsAlpha = -2
    assert excinfo.value.args[0] == "The maximum number of iterations can't be negative"


def test_alpha(vsstls_input_instance):
    vsstls_input_instance.alpha = [-10, 10]
    alpha = vsstls_input_instance.alpha
    assert all(x == y for x, y in zip(alpha, [-10, 10]))
    for a in [[-1.0], [1, 2, 3], [10, -10]]:
        with pytest.raises(RuntimeError) as excinfo:
            vsstls_input_instance.alpha = a
        assert excinfo.value.args[0] == "Invalid guess for free parameter calculation"


def test_couplingResolution(vsstls_input_instance):
    vsstls_input_instance.couplingResolution = 0.01
    couplingResolution = vsstls_input_instance.couplingResolution
    assert couplingResolution == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        vsstls_input_instance.couplingResolution = -0.1
    assert (
        excinfo.value.args[0]
        == "The coupling parameter resolution must be larger than zero"
    )


def test_degeneracyResolution(vsstls_input_instance):
    vsstls_input_instance.degeneracyResolution = 0.01
    degeneracyResolution = vsstls_input_instance.degeneracyResolution
    assert degeneracyResolution == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        vsstls_input_instance.degeneracyResolution = -0.1
    assert (
        excinfo.value.args[0]
        == "The degeneracy parameter resolution must be larger than zero"
    )


def test_freeEnergyIntegrand(vsstls_input_instance):
    arr1 = np.zeros(10)
    arr2 = np.zeros((3, 10))
    fxc = qp.FreeEnergyIntegrand()
    fxc.grid = arr1
    fxc.integrand = arr2
    vsstls_input_instance.freeEnergyIntegrand = fxc
    assert np.array_equal(arr1, vsstls_input_instance.freeEnergyIntegrand.grid)
    assert np.array_equal(arr2, vsstls_input_instance.freeEnergyIntegrand.integrand)
    with pytest.raises(RuntimeError) as excinfo:
        arr1 = np.zeros(10)
        arr2 = np.zeros((2, 10))
        fxc = qp.FreeEnergyIntegrand()
        fxc.grid = arr1
        fxc.integrand = arr2
        vsstls_input_instance.freeEnergyIntegrand = fxc
    assert (
        excinfo.value.args[0]
        == "The free energy integrand does not contain enough temperature points"
    )
    with pytest.raises(RuntimeError) as excinfo:
        arr1 = np.zeros(10)
        arr2 = np.zeros((3, 11))
        fxc = qp.FreeEnergyIntegrand()
        fxc.grid = arr1
        fxc.integrand = arr2
        vsstls_input_instance.freeEnergyIntegrand = fxc
    assert excinfo.value.args[0] == "The free energy integrand is inconsistent"
    with pytest.raises(RuntimeError) as excinfo:
        arr1 = np.zeros(2)
        arr2 = np.zeros((3, 2))
        fxc = qp.FreeEnergyIntegrand()
        fxc.grid = arr1
        fxc.integrand = arr2
        vsstls_input_instance.freeEnergyIntegrand = fxc
    assert (
        excinfo.value.args[0]
        == "The free energy integrand does not contain enough points"
    )


def test_isEqual(vsstls_input_instance):
    thisVSStls = qp.VSStlsInput()
    assert vsstls_input_instance.isEqual(thisVSStls)
    thisVSStls.coupling = 2.0
    thisVSStls.theory = "QSTLS"
    assert not vsstls_input_instance.isEqual(thisVSStls)


def test_print(vsstls_input_instance, capfd):
    vsstls_input_instance.print()
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
    assert "Iet mapping scheme = " in captured
    assert "Maximum number of iterations = 0" in captured
    assert "Minimum error for convergence = 0" in captured
    assert "Mixing parameter = 0" in captured
    assert "Output frequency = 0" in captured
    assert "File with recovery data = " in captured
    assert "Guess for the free parameter = 0,0" in captured
    assert "Resolution for the coupling parameter grid = 0" in captured
    assert "Resolution for the degeneracy parameter grid = 0" in captured
    assert "Minimum error for convergence (alpha) = 0" in captured
    assert "Maximum number of iterations (alpha) = 0" in captured
