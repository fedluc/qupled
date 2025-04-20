import numpy as np
import pytest
from qupled.vsstls import VSStls, Input

import qupled.native as native
import qupled.stls as stls
import qupled.vsstls as vsstls


@pytest.fixture
def scheme():
    scheme = vsstls.VSStls()
    return scheme


def test_vsstls_inheritance():
    assert issubclass(vsstls.VSStls, stls.Stls)


def test_vsstls_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.__init__")
    scheme = vsstls.VSStls()
    super_init.assert_called_once()
    assert isinstance(scheme.results, vsstls.Result)
    assert scheme.native_scheme_cls == native.VSStls
    assert isinstance(scheme.native_inputs, native.VSStlsInput)


def test_compute(mocker, scheme):
    fill_free_energy_integrand = mocker.patch(
        "qupled.vsstls.VSStls._fill_free_energy_integrand"
    )
    super_compute = mocker.patch("qupled.stls.Stls.compute")
    inputs = mocker.ANY
    scheme.compute(inputs)
    fill_free_energy_integrand.assert_called_once_with(inputs)
    super_compute.assert_called_once_with(inputs)


def test_fill_free_energy_integrand(mocker, scheme):
    get_missing_state_points = mocker.patch(
        "qupled.vsstls.VSStls._get_missing_state_points"
    )
    compute = mocker.patch("qupled.vsstls.VSStls.compute")
    update_input_data = mocker.patch("qupled.vsstls.VSStls._update_input_data")
    inputs = mocker.Mock()
    inputs.coupling = mocker.ANY
    inputs.theory = mocker.ANY
    missing_state_points = [0.1, 0.2, 0.3]
    get_missing_state_points.return_value = missing_state_points
    scheme._fill_free_energy_integrand(inputs)
    get_missing_state_points.assert_called_once_with(inputs)
    assert compute.call_count == len(missing_state_points)
    assert update_input_data.call_count == len(missing_state_points)
    for coupling in missing_state_points:
        compute.assert_any_call(inputs)
        update_input_data.assert_any_call(inputs)
    assert inputs.coupling == mocker.ANY


def test_get_missing_state_points_with_no_actual_grid(mocker):
    inputs = mocker.Mock()
    coupling = 1.0
    coupling_resolution = 0.1
    inputs.coupling = coupling
    inputs.coupling_resolution = coupling_resolution
    inputs.free_energy_integrand = mocker.Mock()
    inputs.free_energy_integrand.grid = None
    expected_grid = np.arange(
        coupling_resolution,
        coupling - 0.1 * coupling_resolution,
        3 * coupling_resolution,
    )
    result = VSStls._get_missing_state_points(inputs)
    np.testing.assert_array_equal(result, expected_grid)


def test_get_missing_state_points_with_actual_grid(mocker):
    inputs = mocker.Mock(spec=Input)
    coupling = 1.0
    coupling_resolution = 0.1
    inputs.coupling = coupling
    inputs.coupling_resolution = coupling_resolution
    inputs.free_energy_integrand = mocker.Mock()
    inputs.free_energy_integrand.grid = np.array([0.1, 0.4, 0.7])
    expected_grid = np.arange(
        coupling_resolution,
        coupling - 0.1 * coupling_resolution,
        3 * coupling_resolution,
    )
    precision = int(np.abs(np.log10(0.1)))
    missing_points = np.setdiff1d(
        np.round(expected_grid, precision),
        np.round(inputs.free_energy_integrand.grid, precision),
    )
    result = VSStls._get_missing_state_points(inputs)
    np.testing.assert_array_equal(result, missing_points)


def test_update_input_data(mocker, scheme):
    free_energy_integrand_mock = mocker.patch("qupled.vsstls.FreeEnergyIntegrand")
    inputs = mocker.Mock()
    scheme.results.free_energy_grid = mocker.ANY
    scheme.results.free_energy_integrand = mocker.ANY
    scheme.results.alpha = mocker.ANY
    scheme._update_input_data(inputs)
    free_energy_integrand_mock.assert_called_once_with(
        scheme.results.free_energy_grid,
        scheme.results.free_energy_integrand,
        scheme.results.alpha,
    )
    assert inputs.free_energy_integrand == free_energy_integrand_mock.return_value


def test_get_free_energy_ingtegrand_with_default_database_name(mocker):
    read_results = mocker.patch("qupled.output.DataBase.read_results")
    run_id = mocker.ANY
    read_results.return_value = {
        "free_energy_grid": mocker.ANY,
        "free_energy_integrand": mocker.ANY,
        "alpha": mocker.ANY,
    }
    fxci = vsstls.VSStls.get_free_energy_integrand(run_id)
    assert fxci.grid == read_results.return_value["free_energy_grid"]
    assert fxci.integrand == read_results.return_value["free_energy_integrand"]
    assert fxci.alpha == read_results.return_value["alpha"]
    read_results.assert_called_once_with(
        run_id, None, ["free_energy_grid", "free_energy_integrand", "alpha"]
    )


def test_get_free_energy_ingtegrand_with_custom_database_name(mocker):
    read_results = mocker.patch("qupled.output.DataBase.read_results")
    run_id = mocker.ANY
    database_name = mocker.ANY
    read_results.return_value = {
        "free_energy_grid": mocker.ANY,
        "free_energy_integrand": mocker.ANY,
        "alpha": mocker.ANY,
    }
    fxci = vsstls.VSStls.get_free_energy_integrand(run_id, database_name)
    assert fxci.grid == read_results.return_value["free_energy_grid"]
    assert fxci.integrand == read_results.return_value["free_energy_integrand"]
    assert fxci.alpha == read_results.return_value["alpha"]
    read_results.assert_called_once_with(
        run_id, None, ["free_energy_grid", "free_energy_integrand", "alpha"]
    )


def test_vsstls_input_inheritance():
    assert issubclass(vsstls.Input, stls.Input)


def test_vsstls_input_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Input.__init__")
    free_energy_integrand = mocker.patch("qupled.vsstls.FreeEnergyIntegrand")
    coupling = mocker.ANY
    degeneracy = mocker.ANY
    input = vsstls.Input(coupling, degeneracy)
    assert input.alpha == [0.5, 1.0]
    assert input.coupling_resolution == 0.1
    assert input.degeneracy_resolution == 0.1
    assert input.error_alpha == 1.0e-3
    assert input.iterations_alpha == 50
    assert input.free_energy_integrand == free_energy_integrand.return_value
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "VSSTLS"


def test_vsstls_result_inheritance():
    assert issubclass(vsstls.Result, vsstls.Result)


def test_vsstls_result_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Result.__init__")
    stls_results = vsstls.Result()
    assert stls_results.free_energy_grid is None
    assert stls_results.free_energy_integrand is None
    assert stls_results.alpha is None
    super_init.assert_called_once()


def test_free_energy_integrand_initialization(mocker):
    grid = mocker.ANY
    integrand = mocker.ANY
    fxci = vsstls.FreeEnergyIntegrand(grid, integrand)
    assert fxci.grid == grid
    assert fxci.integrand == integrand


def test_free_energy_integrand_initialization_defaults():
    fxci = vsstls.FreeEnergyIntegrand()
    assert fxci.grid is None
    assert fxci.integrand is None


def test_free_energy_integrand_to_native(mocker):
    FreeEnergyIntegrand = mocker.patch("qupled.native.FreeEnergyIntegrand")
    native_fxci = mocker.ANY
    grid = mocker.ANY
    integrand = mocker.ANY
    FreeEnergyIntegrand.return_value = native_fxci
    guess = vsstls.FreeEnergyIntegrand(grid, integrand)
    result = guess.to_native()
    assert result == native_fxci
    assert result.grid == grid
    assert result.integrand == integrand
