import pytest
from qupled.vsstls import VSStls
from qupled.stls import Stls
from qupled.native import VSStls as NativeVSStls
from qupled.base import IterativeScheme, Input, Result


def test_vsstls_inheritance():
    assert issubclass(VSStls, IterativeScheme)


def test_vsstls_compute(mocker):
    super_compute = mocker.patch("qupled.base.ClassicScheme.compute")
    input = mocker.ANY
    scheme = VSStls()
    scheme.compute(input)
    super_compute.assert_called_once_with(input, NativeVSStls, mocker.ANY, mocker.ANY)


def test_get_free_energy_integrand(mocker):
    run_id = mocker.ANY
    database_name = mocker.ANY
    names = mocker.ANY
    mock_read_results = mocker.patch("qupled.util.DataBase.read_results")
    mock_read_results.return_value = {
        "free_energy_grid": [1, 2, 3],
        "free_energy_integrand": [0.1, 0.2, 0.3],
        "alpha": [0.5, 1.0],
    }
    mock_free_energy_integrand = mocker.patch("qupled.native.FreeEnergyIntegrand")
    mock_fxci_instance = mock_free_energy_integrand.return_value
    result = VSStls.get_free_energy_integrand(run_id, database_name)
    mock_read_results.assert_called_once_with(run_id, database_name, names)
    mock_free_energy_integrand.assert_called_once()
    assert result == mock_fxci_instance
    assert result.grid == [1, 2, 3]
    assert result.integrand == [0.1, 0.2, 0.3]
    assert result.alpha == [0.5, 1.0]


def test_vsstls_input_inheritance():
    assert issubclass(VSStls.Input, Stls.Input)


def test_vsstls_input_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.Input.__init__")
    free_energy_integrand = mocker.patch("qupled.native.FreeEnergyIntegrand")
    coupling = mocker.ANY
    degeneracy = mocker.ANY
    input = VSStls.Input(coupling, degeneracy)
    assert input.alpha == [0.5, 1.0]
    assert input.coupling_resolution == 0.1
    assert input.degeneracy_resolution == 0.1
    assert input.error_alpha == 1.0e-3
    assert input.iterations_alpha == 50
    assert input.free_energy_integrand == free_energy_integrand.return_value
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "VSSTLS"


def test_vsstls_result_inheritance():
    assert issubclass(VSStls.Result, Stls.Result)


def test_vsstls_result_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.Result.__init__")
    stls_results = VSStls.Result()
    assert stls_results.free_energy_grid is None
    assert stls_results.free_energy_integrand is None
    assert stls_results.alpha is None
    super_init.assert_called_once()
