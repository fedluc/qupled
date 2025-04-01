import qupled.native as native
import qupled.stls as stls
import qupled.vsstls as vsstls


def test_vsstls_inheritance():
    assert issubclass(vsstls.VSStls, stls.Stls)


def test_stls_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.__init__")
    scheme = vsstls.VSStls()
    super_init.assert_called_once()
    assert isinstance(scheme.results, vsstls.Result)
    assert scheme.native_scheme_cls == native.VSStls
    assert isinstance(scheme.native_inputs, native.VSStlsInput)


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
    alpha = mocker.ANY
    fxci = vsstls.FreeEnergyIntegrand(grid, integrand, alpha)
    assert fxci.grid == grid
    assert fxci.integrand == integrand
    assert fxci.alpha == alpha


def test_free_energy_integrand_initialization_defaults():
    fxci = vsstls.FreeEnergyIntegrand()
    assert fxci.grid is None
    assert fxci.integrand is None
    assert fxci.alpha is None


def test_free_energy_integrand_to_native(mocker):
    FreeEnergyIntegrand = mocker.patch("qupled.native.FreeEnergyIntegrand")
    native_fxci = mocker.ANY
    grid = mocker.ANY
    integrand = mocker.ANY
    alpha = mocker.ANY
    FreeEnergyIntegrand.return_value = native_fxci
    guess = vsstls.FreeEnergyIntegrand(grid, integrand, alpha)
    result = guess.to_native()
    assert result == native_fxci
    assert result.grid == grid
    assert result.integrand == integrand
    assert result.alpha == alpha
