import numpy as np
import qupled.native as native
import qupled.stls as stls
import qupled.rpa as rpa


def test_stls_inheritance():
    assert issubclass(stls.Stls, rpa.Rpa)


def test_stls_initialization(mocker):
    super_init = mocker.patch("qupled.rpa.Rpa.__init__")
    scheme = stls.Stls()
    super_init.assert_called_once()
    assert isinstance(scheme.results, stls.Result)
    assert scheme.native_scheme_cls == native.Stls
    assert isinstance(scheme.native_inputs, native.StlsInput)


def test_get_initial_guess_with_default_database_name(mocker):
    read_results = mocker.patch("qupled.output.DataBase.read_results")
    run_id = mocker.ANY
    read_results.return_value = {
        "wvg": np.array([1.0, 2.0, 3.0]),
        "slfc": np.array([0.1, 0.2, 0.3]),
    }
    guess = stls.Stls.get_initial_guess(run_id)
    assert np.array_equal(guess.wvg, read_results.return_value["wvg"])
    assert np.array_equal(guess.slfc, read_results.return_value["slfc"])
    read_results.assert_called_once_with(run_id, None, ["wvg", "slfc"])


def test_get_initial_guess_with_custom_database_name(mocker):
    read_results = mocker.patch("qupled.output.DataBase.read_results")
    run_id = mocker.ANY
    database_name = mocker.ANY
    read_results.return_value = {
        "wvg": np.array([1.0, 2.0, 3.0]),
        "slfc": np.array([0.1, 0.2, 0.3]),
    }
    guess = stls.Stls.get_initial_guess(run_id, database_name)
    assert np.array_equal(guess.wvg, read_results.return_value["wvg"])
    assert np.array_equal(guess.slfc, read_results.return_value["slfc"])
    read_results.assert_called_once_with(run_id, database_name, ["wvg", "slfc"])


def test_stls_input_inheritance():
    assert issubclass(stls.Input, rpa.Input)


def test_stls_input_initialization(mocker):
    super_init = mocker.patch("qupled.rpa.Input.__init__")
    guess = mocker.patch("qupled.stls.Guess")
    coupling = 1.5
    degeneracy = 3.0
    input = stls.Input(coupling, degeneracy)
    assert input.error == 1.0e-5
    assert input.mixing == 1.0
    assert input.iterations == 1000
    assert input.output_frequency == 10
    assert input.guess == guess.return_value
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "STLS"


def test_stls_result_inheritance():
    assert issubclass(stls.Result, rpa.Result)


def test_stls_result_initialization(mocker):
    super_init = mocker.patch("qupled.rpa.Result.__init__")
    results = stls.Result()
    assert results.error is None
    super_init.assert_called_once()


def test_stls_guess_initialization(mocker):
    wvg = mocker.ANY
    slfc = mocker.ANY
    guess = stls.Guess(wvg, slfc)
    assert guess.wvg == wvg
    assert guess.slfc == slfc


def test_stls_guess_initialization_defaults():
    guess = stls.Guess()
    assert guess.wvg is None
    assert guess.slfc is None


def test_stls_guess_to_native(mocker):
    StlsGuess = mocker.patch("qupled.native.StlsGuess")
    native_guess = mocker.ANY
    wvg = mocker.ANY
    slfc = mocker.ANY
    StlsGuess.return_value = native_guess
    guess = stls.Guess(wvg, slfc)
    result = guess.to_native()
    assert result == native_guess
    assert result.wvg == wvg
    assert result.slfc == slfc
