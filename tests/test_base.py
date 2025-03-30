import numpy as np
from qupled.base import (
    Input,
    Result,
    ClassicScheme,
    IterativeScheme,
    QuantumIterativeScheme,
)
from qupled.database import DataBaseHandler
from unittest import main
import pytest


class TestQuantumIterativeSchemeGetInitialGuess:

    def test_get_initial_guess_with_default_database_name(self, mocker):
        read_run = mocker.patch("qupled.util.DataBase.read_run")
        run_id = "test_run"
        read_run.return_value = {
            DataBaseHandler.INPUTS_TABLE_NAME: {"matsubara": 128},
            DataBaseHandler.RESULTS_TABLE_NAME: {
                "wvg": np.array([1, 2, 3]),
                "ssf": np.array([4, 5, 6]),
                "adr": np.array([7, 8, 9]),
            },
        }
        guess = QuantumIterativeScheme.get_initial_guess(run_id)
        assert (guess.wvg == np.array([1, 2, 3])).all()
        assert (guess.ssf == np.array([4, 5, 6])).all()
        assert (guess.adr == np.array([7, 8, 9])).all()
        assert guess.matsubara == 128
        read_run.assert_called_once_with(
            run_id, None, ["matsubara"], ["wvg", "ssf", "adr"]
        )

    def test_get_initial_guess_with_custom_database_name(self, mocker):
        read_run = mocker.patch("qupled.util.DataBase.read_run")
        run_id = "test_run"
        database_name = "test_db"
        read_run.return_value = {
            DataBaseHandler.INPUTS_TABLE_NAME: {"matsubara": 128},
            DataBaseHandler.RESULTS_TABLE_NAME: {
                "wvg": np.array([1, 2, 3]),
                "ssf": np.array([4, 5, 6]),
                "adr": np.array([7, 8, 9]),
            },
        }
        guess = QuantumIterativeScheme.get_initial_guess(run_id, database_name)
        assert (guess.wvg == np.array([1, 2, 3])).all()
        assert (guess.ssf == np.array([4, 5, 6])).all()
        assert (guess.adr == np.array([7, 8, 9])).all()
        assert guess.matsubara == 128
        read_run.assert_called_once_with(
            run_id, database_name, ["matsubara"], ["wvg", "ssf", "adr"]
        )


class TestQuantumIterativeSchemeGuess:

    def test_initialization(self):
        guess = QuantumIterativeScheme.Guess(
            wvg=np.array([1, 2, 3]),
            ssf=np.array([4, 5, 6]),
            adr=np.array([7, 8, 9]),
            matsubara=128,
        )
        assert (guess.wvg == np.array([1, 2, 3])).all()
        assert (guess.ssf == np.array([4, 5, 6])).all()
        assert (guess.adr == np.array([7, 8, 9])).all()
        assert guess.matsubara == 128

    def test_initialization_defaults(self):
        guess = QuantumIterativeScheme.Guess()
        assert guess.wvg is None
        assert guess.ssf is None
        assert guess.adr is None
        assert guess.matsubara == 0

    def test_to_native(self, mocker):
        QstlsGuess = mocker.patch("qupled.base.native.QstlsGuess")
        native_guess = mocker.ANY
        QstlsGuess.return_value = native_guess
        guess = QuantumIterativeScheme.Guess(
            wvg=np.array([1, 2, 3]),
            ssf=np.array([4, 5, 6]),
            adr=np.array([7, 8, 9]),
            matsubara=128,
        )
        result = guess.to_native()
        assert result == native_guess
        assert (result.wvg == np.array([1, 2, 3])).all()
        assert (result.ssf == np.array([4, 5, 6])).all()
        assert (result.adr == np.array([7, 8, 9])).all()
        assert result.matsubara == 128


if __name__ == "__main__":
    main()
