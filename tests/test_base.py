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


@pytest.fixture
def input_obj():
    return Input(coupling=1.0, degeneracy=2.0)


def test_input_initialization(input_obj):
    assert input_obj.coupling == 1.0
    assert input_obj.degeneracy == 2.0
    assert input_obj.chemical_potential == [-10.0, 10.0]
    assert input_obj.cutoff == 10.0
    assert input_obj.frequency_cutoff == 10.0
    assert input_obj.integral_error == 1.0e-5
    assert input_obj.integral_strategy == "full"
    assert input_obj.matsubara == 128
    assert input_obj.resolution == 0.1
    assert input_obj.threads == 1


def test_input_to_native(mocker, input_obj):
    native_input = mocker.MagicMock()
    input_obj.to_native(native_input)
    assert native_input.coupling == 1.0
    assert native_input.degeneracy == 2.0


@pytest.fixture
def result_obj():
    return Result()


def test_result_initialization(result_obj):
    assert result_obj.idr is None
    assert result_obj.rdf is None
    assert result_obj.rdf_grid is None
    assert result_obj.sdr is None
    assert result_obj.slfc is None
    assert result_obj.ssf is None
    assert result_obj.uint is None
    assert result_obj.wvg is None


def test_result_from_native(mocker, result_obj):
    native_scheme = mocker.MagicMock()
    native_scheme.idr = np.array([1, 2, 3])
    result_obj.from_native(native_scheme)
    assert np.array_equal(result_obj.idr, np.array([1, 2, 3]))


def test_result_compute_rdf_with_default_grid(mocker, result_obj):
    compute_rdf = mocker.patch("qupled.base.native.compute_rdf")
    result_obj.wvg = np.array([1, 2, 3])
    result_obj.ssf = np.array([4, 5, 6])
    compute_rdf.return_value = np.array([7, 8, 9])
    result_obj.compute_rdf()
    assert result_obj.rdf is not None
    assert result_obj.rdf_grid is not None
    compute_rdf.assert_called_once()
    assert np.array_equal(result_obj.rdf, np.array([7, 8, 9]))
    assert np.array_equal(result_obj.rdf_grid, np.arange(0.0, 10.0, 0.01))


def test_result_compute_rdf_with_custom_grid(mocker, result_obj):
    compute_rdf = mocker.patch("qupled.base.native.compute_rdf")
    rdf_grid = np.array([0, 1, 2])
    result_obj.wvg = np.array([1, 2, 3])
    result_obj.ssf = np.array([4, 5, 6])
    compute_rdf.return_value = np.array([7, 8, 9])
    result_obj.compute_rdf(rdf_grid)
    assert result_obj.rdf is not None
    assert result_obj.rdf_grid is not None
    compute_rdf.assert_called_once()
    assert np.array_equal(result_obj.rdf, np.array([7, 8, 9]))
    assert np.array_equal(result_obj.rdf_grid, rdf_grid)


@pytest.fixture
def classic_scheme(mocker):
    scheme = ClassicScheme()
    scheme.db_handler = mocker.MagicMock()
    return scheme


class TestClassicScheme:

    def test_run_id(self, classic_scheme):
        run_id = "run_id"
        classic_scheme.db_handler.run_id = run_id
        assert classic_scheme.run_id == run_id

    def test_compute(self, classic_scheme, mocker):
        inputs = Input(coupling=1.0, degeneracy=2.0)
        native_input = mocker.MagicMock()
        native_scheme = mocker.MagicMock()
        native_scheme.compute.return_value = 0
        native_scheme.recovery = "recovery_data"
        native_scheme_cls = mocker.MagicMock(return_value=native_scheme)
        native_scheme_cls.compute.return_value = 0
        result = Result()
        classic_scheme.compute(inputs, native_scheme_cls, native_input, result)
        assert classic_scheme.inputs is not None
        assert classic_scheme.results is not None
        native_scheme_cls.assert_called_once_with(native_input)
        native_scheme.compute.assert_called_once()
        classic_scheme.db_handler.insert_run.assert_called_once_with(inputs, result)

    def test_check_status_and_clean_success(self, classic_scheme, mocker):
        exists = mocker.patch("os.path.exists")
        remove = mocker.patch("os.remove")
        exit = mocker.patch("sys.exit")
        exists.return_value = True
        classic_scheme._check_status_and_clean(0, "recovery_file")
        exists.assert_called_once_with("recovery_file")
        remove.assert_called_once_with("recovery_file")
        exit.assert_not_called()

    def test_check_status_and_clean_failure(self, classic_scheme, mocker):
        exists = mocker.patch("os.path.exists")
        exit = mocker.patch("sys.exit")
        exists.return_value = False
        classic_scheme._check_status_and_clean(1, "recovery_file")
        exists.assert_not_called()
        exit.assert_called_once_with("Error while solving the dielectric theory")

    def test_save(self, classic_scheme):
        classic_scheme.results = Result()
        classic_scheme._save()
        classic_scheme.db_handler.insert_run.assert_called_once()

    def test_save_with_results(self, classic_scheme):
        classic_scheme.results = Result()
        classic_scheme.inputs = Input(coupling=1.0, degeneracy=2.0)
        classic_scheme._save()
        classic_scheme.db_handler.insert_run.assert_called_once_with(
            classic_scheme.inputs, classic_scheme.results
        )

    def test_save_without_results(self, classic_scheme):
        classic_scheme.results = None
        classic_scheme._save()
        classic_scheme.db_handler.insert_run.assert_not_called()

    def test_compute_rdf_with_default_grid(self, classic_scheme, mocker):
        compute_rdf = mocker.patch("qupled.base.Result.compute_rdf")
        classic_scheme.results = Result()
        classic_scheme.compute_rdf()
        compute_rdf.assert_called_once_with(None)
        classic_scheme.db_handler.insert_results.assert_called_once_with(
            {
                "rdf": classic_scheme.results.rdf,
                "rdf_grid": classic_scheme.results.rdf_grid,
            }
        )

    def test_compute_rdf_with_custom_grid(self, classic_scheme, mocker):
        compute_rdf = mocker.patch("qupled.base.Result.compute_rdf")
        classic_scheme.results = Result()
        rdf_grid = np.array([1, 2, 3])
        classic_scheme.compute_rdf(rdf_grid)
        compute_rdf.assert_called_once_with(rdf_grid)
        classic_scheme.db_handler.insert_results.assert_called_once_with(
            {
                "rdf": classic_scheme.results.rdf,
                "rdf_grid": classic_scheme.results.rdf_grid,
            }
        )

    def test_compute_rdf_without_results(self, classic_scheme):
        classic_scheme.results = None
        classic_scheme.compute_rdf()
        classic_scheme.db_handler.insert_results.assert_not_called()


class TestIterativeScheme:

    def test_get_initial_guess_with_default_database_name(self, mocker):
        read_results = mocker.patch("qupled.util.DataBase.read_results")
        run_id = "test_run"
        read_results.return_value = {
            "wvg": np.array([1, 2, 3]),
            "slfc": np.array([4, 5, 6]),
        }
        guess = IterativeScheme.get_initial_guess(run_id)
        assert (guess.wvg == np.array([1, 2, 3])).all()
        assert (guess.slfc == np.array([4, 5, 6])).all()
        read_results.assert_called_once_with(run_id, None, ["wvg", "slfc"])

    def test_get_initial_guess_with_custom_database_name(self, classic_scheme, mocker):
        read_results = mocker.patch("qupled.util.DataBase.read_results")
        run_id = "test_run"
        database_name = "test_db"
        read_results.return_value = {
            "wvg": np.array([1, 2, 3]),
            "slfc": np.array([4, 5, 6]),
        }
        guess = IterativeScheme.get_initial_guess(run_id, database_name)
        assert (guess.wvg == np.array([1, 2, 3])).all()
        assert (guess.slfc == np.array([4, 5, 6])).all()
        read_results.assert_called_once_with(run_id, database_name, ["wvg", "slfc"])


class TestIterativeSchemeGuess:

    def test_initialization(self):
        guess = IterativeScheme.Guess(wvg=np.array([1, 2, 3]), slfc=np.array([4, 5, 6]))
        assert (guess.wvg == np.array([1, 2, 3])).all()
        assert (guess.slfc == np.array([4, 5, 6])).all()

    def test_initialization_defaults(self):
        guess = IterativeScheme.Guess()
        assert guess.wvg is None
        assert guess.slfc is None

    def test_to_native(self, mocker):
        StlsGuess = mocker.patch("qupled.base.native.StlsGuess")
        native_guess = mocker.MagicMock()
        StlsGuess.return_value = native_guess
        guess = IterativeScheme.Guess(wvg=np.array([1, 2, 3]), slfc=np.array([4, 5, 6]))
        result = guess.to_native()
        assert result == native_guess
        assert (result.wvg == np.array([1, 2, 3])).all()
        assert (result.slfc == np.array([4, 5, 6])).all()


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
        native_guess = mocker.MagicMock()
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
