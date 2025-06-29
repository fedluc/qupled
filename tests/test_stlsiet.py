import numpy as np
import pytest

import qupled.stls as stls
import qupled.stlsiet as stlsiet
from qupled.database import DataBaseHandler


def test_stls_iet_inheritance():
    assert issubclass(stlsiet.Solver, stls.Solver)


def test_stls_iet_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Solver.__init__")
    scheme = stlsiet.Solver()
    super_init.assert_called_once()
    assert isinstance(scheme.results, stlsiet.Result)


def test_stls_iet_input_inheritance():
    assert issubclass(stlsiet.Input, stls.Input)


def test_stls_iet_input_initialization_valid_theory(mocker):
    theory = "STLS-HNC"
    input = stlsiet.Input(mocker.ANY, mocker.ANY, theory=theory)
    assert input.theory == theory


def test_stls_iet_input_initialization_invalid_theory():
    with pytest.raises(ValueError):
        stlsiet.Input(1.0, 1.0, "INVALID-THEORY")


def test_stls_iet_input_initialization_default_theory():
    with pytest.raises(ValueError):
        stlsiet.Input(1.0, 1.0)


def test_stls_iet_result_inheritance():
    assert issubclass(stlsiet.Result, stls.Result)


def test_stls_iet_result_initialization(mocker):
    results = stlsiet.Result()
    assert results.bf is None


def test_get_initial_guess_with_default_database_name(mocker):
    read_results = mocker.patch("qupled.output.DataBase.read_results")
    run_id = mocker.ANY
    read_results.return_value = {
        "wvg": np.array([1, 2, 3]),
        "ssf": np.array([4, 5, 6]),
        "lfc": np.array([7, 8, 9]),
    }
    guess = stlsiet.Solver.get_initial_guess(run_id)
    assert np.array_equal(guess.wvg, np.array([1, 2, 3]))
    assert np.array_equal(guess.ssf, np.array([4, 5, 6]))
    assert np.array_equal(guess.lfc, np.array([7, 8, 9]))
    read_results.assert_called_once_with(run_id, None, ["wvg", "ssf", "lfc"])


def test_get_initial_guess_with_custom_database_name(mocker):
    read_results = mocker.patch("qupled.output.DataBase.read_results")
    database_name = mocker.ANY
    run_id = mocker.ANY
    read_results.return_value = {
        "wvg": np.array([1, 2, 3]),
        "ssf": np.array([4, 5, 6]),
        "lfc": np.array([7, 8, 9]),
    }
    guess = stlsiet.Solver.get_initial_guess(run_id)
    assert np.array_equal(guess.wvg, np.array([1, 2, 3]))
    assert np.array_equal(guess.ssf, np.array([4, 5, 6]))
    assert np.array_equal(guess.lfc, np.array([7, 8, 9]))
    read_results.assert_called_once_with(run_id, database_name, ["wvg", "ssf", "lfc"])
