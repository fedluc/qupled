import pytest
import numpy as np
from qupled.native import QVSStlsInput, VSInput, FreeEnergyIntegrand


@pytest.fixture
def qvsstls_input_instance():
    return QVSStlsInput()


def test_init(qvsstls_input_instance):
    assert issubclass(QVSStlsInput, VSInput)
    assert hasattr(qvsstls_input_instance, "error_alpha")
    assert hasattr(qvsstls_input_instance, "iterations_alpha")
    assert hasattr(qvsstls_input_instance, "alpha")
    assert hasattr(qvsstls_input_instance, "coupling_resolution")
    assert hasattr(qvsstls_input_instance, "degeneracy_resolution")
    assert hasattr(qvsstls_input_instance, "free_energy_integrand")
    assert hasattr(qvsstls_input_instance, "guess")
    assert hasattr(qvsstls_input_instance.guess, "wvg")
    assert hasattr(qvsstls_input_instance.guess, "ssf")
    assert hasattr(qvsstls_input_instance.guess, "adr")
    assert hasattr(qvsstls_input_instance.guess, "matsubara")
    assert hasattr(qvsstls_input_instance, "fixed")
    assert hasattr(qvsstls_input_instance.free_energy_integrand, "grid")
    assert hasattr(qvsstls_input_instance.free_energy_integrand, "integrand")


def test_defaults(qvsstls_input_instance):
    assert np.isnan(qvsstls_input_instance.error_alpha)
    assert qvsstls_input_instance.iterations_alpha == 0
    assert qvsstls_input_instance.alpha.size == 0
    assert np.isnan(qvsstls_input_instance.coupling_resolution)
    assert np.isnan(qvsstls_input_instance.degeneracy_resolution)
    assert qvsstls_input_instance.free_energy_integrand.grid.size == 0
    assert qvsstls_input_instance.free_energy_integrand.integrand.size == 0
    assert qvsstls_input_instance.guess.wvg.size == 0
    assert qvsstls_input_instance.guess.ssf.size == 0
    assert qvsstls_input_instance.guess.adr.size == 0
    assert qvsstls_input_instance.guess.matsubara == 0
    assert qvsstls_input_instance.fixed == ""


def test_fixed(qvsstls_input_instance):
    qvsstls_input_instance.fixed = "fixedFile"
    fixed = qvsstls_input_instance.fixed
    assert fixed == "fixedFile"


def test_error_alpha(qvsstls_input_instance):
    qvsstls_input_instance.error_alpha = 0.001
    error_alpha = qvsstls_input_instance.error_alpha
    assert error_alpha == 0.001
    with pytest.raises(RuntimeError) as excinfo:
        qvsstls_input_instance.error_alpha = -0.1
    assert (
        excinfo.value.args[0]
        == "The minimum error for convergence must be larger than zero"
    )


def test_iterations_alpha(qvsstls_input_instance):
    qvsstls_input_instance.iterations_alpha = 1
    iterations_alpha = qvsstls_input_instance.iterations_alpha
    assert iterations_alpha == 1
    with pytest.raises(RuntimeError) as excinfo:
        qvsstls_input_instance.iterations_alpha = -2
    assert excinfo.value.args[0] == "The maximum number of iterations can't be negative"


def test_alpha(qvsstls_input_instance):
    qvsstls_input_instance.alpha = [-10, 10]
    alpha = qvsstls_input_instance.alpha
    assert all(x == y for x, y in zip(alpha, [-10, 10]))
    for a in [[-1.0], [1, 2, 3], [10, -10]]:
        with pytest.raises(RuntimeError) as excinfo:
            qvsstls_input_instance.alpha = a
        assert excinfo.value.args[0] == "Invalid guess for free parameter calculation"


def test_coupling_resolution(qvsstls_input_instance):
    qvsstls_input_instance.coupling_resolution = 0.01
    coupling_resolution = qvsstls_input_instance.coupling_resolution
    assert coupling_resolution == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        qvsstls_input_instance.coupling_resolution = -0.1
    assert (
        excinfo.value.args[0]
        == "The coupling parameter resolution must be larger than zero"
    )


def test_degeneracy_resolution(qvsstls_input_instance):
    qvsstls_input_instance.degeneracy_resolution = 0.01
    degeneracy_resolution = qvsstls_input_instance.degeneracy_resolution
    assert degeneracy_resolution == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        qvsstls_input_instance.degeneracy_resolution = -0.1
    assert (
        excinfo.value.args[0]
        == "The degeneracy parameter resolution must be larger than zero"
    )


def test_free_energy_integrand(qvsstls_input_instance):
    arr1 = np.zeros(10)
    arr2 = np.zeros((3, 10))
    fxc = FreeEnergyIntegrand()
    fxc.grid = arr1
    fxc.integrand = arr2
    qvsstls_input_instance.free_energy_integrand = fxc
    assert np.array_equal(arr1, qvsstls_input_instance.free_energy_integrand.grid)
    assert np.array_equal(arr2, qvsstls_input_instance.free_energy_integrand.integrand)


def test_free_energy_integrand_inconsistent(qvsstls_input_instance):
    with pytest.raises(RuntimeError) as excinfo:
        arr1 = np.zeros(10)
        arr2 = np.zeros((3, 11))
        fxc = FreeEnergyIntegrand()
        fxc.grid = arr1
        fxc.integrand = arr2
        qvsstls_input_instance.free_energy_integrand = fxc
    assert excinfo.value.args[0] == "The free energy integrand is inconsistent"


def test_is_equal_default(qvsstls_input_instance):
    assert not qvsstls_input_instance.is_equal(qvsstls_input_instance)


def test_is_equal(qvsstls_input_instance):
    this_qvsstls = QVSStlsInput()
    this_qvsstls.coupling = 2.0
    this_qvsstls.degeneracy = 1.0
    this_qvsstls.integral_error = 0.1
    this_qvsstls.threads = 1
    this_qvsstls.theory = "STLS"
    this_qvsstls.matsubara = 1
    this_qvsstls.resolution = 0.1
    this_qvsstls.cutoff = 1.0
    this_qvsstls.error = 0.1
    this_qvsstls.mixing = 1.0
    this_qvsstls.output_frequency = 1
    this_qvsstls.coupling_resolution = 0.1
    this_qvsstls.degeneracy_resolution = 0.1
    this_qvsstls.error_alpha = 0.1
    this_qvsstls.iterations_alpha = 1
    assert this_qvsstls.is_equal(this_qvsstls)


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
    assert "Frequency cutoff = nan" in captured
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
