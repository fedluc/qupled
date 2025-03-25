import pytest
import numpy as np
from qupled.native import VSInput, VSStlsInput, FreeEnergyIntegrand


@pytest.fixture
def vsstls_input_instance():
    return VSStlsInput()


def test_init(vsstls_input_instance):
    assert issubclass(VSStlsInput, VSInput)
    assert hasattr(vsstls_input_instance, "error_alpha")
    assert hasattr(vsstls_input_instance, "iterations_alpha")
    assert hasattr(vsstls_input_instance, "alpha")
    assert hasattr(vsstls_input_instance, "coupling_resolution")
    assert hasattr(vsstls_input_instance, "degeneracy_resolution")
    assert hasattr(vsstls_input_instance, "free_energy_integrand")
    assert hasattr(vsstls_input_instance.free_energy_integrand, "grid")
    assert hasattr(vsstls_input_instance.free_energy_integrand, "integrand")
    assert hasattr(vsstls_input_instance, "mapping")
    assert hasattr(vsstls_input_instance, "guess")


def test_defaults(vsstls_input_instance):
    assert np.isnan(vsstls_input_instance.error_alpha)
    assert vsstls_input_instance.iterations_alpha == 0
    assert vsstls_input_instance.alpha.size == 0
    assert np.isnan(vsstls_input_instance.coupling_resolution)
    assert np.isnan(vsstls_input_instance.degeneracy_resolution)
    assert vsstls_input_instance.free_energy_integrand.grid.size == 0
    assert vsstls_input_instance.free_energy_integrand.integrand.size == 0
    assert vsstls_input_instance.guess.wvg.size == 0
    assert vsstls_input_instance.guess.slfc.size == 0


def test_error_alpha(vsstls_input_instance):
    vsstls_input_instance.error_alpha = 0.001
    error_alpha = vsstls_input_instance.error_alpha
    assert error_alpha == 0.001
    with pytest.raises(RuntimeError) as excinfo:
        vsstls_input_instance.error_alpha = -0.1
    assert (
        excinfo.value.args[0]
        == "The minimum error for convergence must be larger than zero"
    )


def test_iterations_alpha(vsstls_input_instance):
    vsstls_input_instance.iterations_alpha = 1
    iterations_alpha = vsstls_input_instance.iterations_alpha
    assert iterations_alpha == 1
    with pytest.raises(RuntimeError) as excinfo:
        vsstls_input_instance.iterations_alpha = -2
    assert excinfo.value.args[0] == "The maximum number of iterations can't be negative"


def test_alpha(vsstls_input_instance):
    vsstls_input_instance.alpha = [-10, 10]
    alpha = vsstls_input_instance.alpha
    assert all(x == y for x, y in zip(alpha, [-10, 10]))
    for a in [[-1.0], [1, 2, 3], [10, -10]]:
        with pytest.raises(RuntimeError) as excinfo:
            vsstls_input_instance.alpha = a
        assert excinfo.value.args[0] == "Invalid guess for free parameter calculation"


def test_coupling_resolution(vsstls_input_instance):
    vsstls_input_instance.coupling_resolution = 0.01
    coupling_resolution = vsstls_input_instance.coupling_resolution
    assert coupling_resolution == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        vsstls_input_instance.coupling_resolution = -0.1
    assert (
        excinfo.value.args[0]
        == "The coupling parameter resolution must be larger than zero"
    )


def test_degeneracy_resolution(vsstls_input_instance):
    vsstls_input_instance.degeneracy_resolution = 0.01
    degeneracy_resolution = vsstls_input_instance.degeneracy_resolution
    assert degeneracy_resolution == 0.01
    with pytest.raises(RuntimeError) as excinfo:
        vsstls_input_instance.degeneracy_resolution = -0.1
    assert (
        excinfo.value.args[0]
        == "The degeneracy parameter resolution must be larger than zero"
    )


def test_free_energy_integrand(vsstls_input_instance):
    arr1 = np.zeros(10)
    arr2 = np.zeros((3, 10))
    fxc = FreeEnergyIntegrand()
    fxc.grid = arr1
    fxc.integrand = arr2
    vsstls_input_instance.free_energy_integrand = fxc
    assert np.array_equal(arr1, vsstls_input_instance.free_energy_integrand.grid)
    assert np.array_equal(arr2, vsstls_input_instance.free_energy_integrand.integrand)


def test_free_energy_integrand_inconsistent(vsstls_input_instance):
    with pytest.raises(RuntimeError) as excinfo:
        arr1 = np.zeros(10)
        arr2 = np.zeros((3, 11))
        fxc = FreeEnergyIntegrand()
        fxc.grid = arr1
        fxc.integrand = arr2
        vsstls_input_instance.free_energy_integrand = fxc
    assert excinfo.value.args[0] == "The free energy integrand is inconsistent"


def test_is_equal_default(vsstls_input_instance):
    assert not vsstls_input_instance.is_equal(vsstls_input_instance)


def test_is_equal(vsstls_input_instance):
    this_vs_stls = VSStlsInput()
    this_vs_stls.coupling = 2.0
    this_vs_stls.degeneracy = 1.0
    this_vs_stls.integral_error = 0.1
    this_vs_stls.threads = 1
    this_vs_stls.theory = "STLS"
    this_vs_stls.matsubara = 1
    this_vs_stls.resolution = 0.1
    this_vs_stls.cutoff = 1.0
    this_vs_stls.error = 0.1
    this_vs_stls.mixing = 1.0
    this_vs_stls.output_frequency = 1
    this_vs_stls.coupling_resolution = 0.1
    this_vs_stls.degeneracy_resolution = 0.1
    this_vs_stls.error_alpha = 0.1
    this_vs_stls.iterations_alpha = 1
    assert this_vs_stls.is_equal(this_vs_stls)


def test_print(vsstls_input_instance, capfd):
    vsstls_input_instance.print()
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
