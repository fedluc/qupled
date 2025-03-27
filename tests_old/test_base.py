import numpy as np
from qupled.base import (
    Input,
    Result,
    ClassicScheme,
    IterativeScheme,
    QuantumIterativeScheme,
)
from qupled.database import DataBaseHandler
from unittest import mock, TestCase, main as unittestmain
import pytest


@pytest.fixture
def input_obj():
    return Input(coupling=1.0, degeneracy=2.0)

def test_initialization(input_obj):
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

def test_to_native(input_obj):
    native_input = mock.MagicMock()
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

def test_result_from_native(result_obj):
    native_scheme = mock.MagicMock()
    native_scheme.idr = np.array([1, 2, 3])
    result_obj.from_native(native_scheme)
    assert np.array_equal(result_obj.idr, np.array([1, 2, 3]))

@mock.patch("qupled.base.native.compute_rdf")
def test_result_compute_rdf_with_default_grid(compute_rdf, result_obj):
    result_obj.wvg = np.array([1, 2, 3])
    result_obj.ssf = np.array([4, 5, 6])
    compute_rdf.return_value = np.array([7, 8, 9])
    result_obj.compute_rdf()
    assert result_obj.rdf is not None
    assert result_obj.rdf_grid is not None
    compute_rdf.assert_called_once()
    assert np.array_equal(result_obj.rdf, np.array([7, 8, 9]))
    assert np.array_equal(result_obj.rdf_grid, np.arange(0.0, 10.0, 0.01))

@mock.patch("qupled.base.native.compute_rdf")
def test_result_compute_rdf_with_custom_grid(compute_rdf, result_obj):
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


class TestClassicScheme(TestCase):

    def setUp(self):
        self.scheme = ClassicScheme()
        self.scheme.db_handler = mock.MagicMock()

    def test_run_id(self):
        self.scheme.db_handler.run_id = "test_id"
        self.assertEqual(self.scheme.run_id, "test_id")

    def test_compute(self):
        inputs = Input(coupling=1.0, degeneracy=2.0)
        native_input = mock.MagicMock()
        native_scheme = mock.MagicMock()
        native_scheme.compute.return_value = 0
        native_scheme.recovery = "recovery_data"
        native_scheme_cls = mock.MagicMock(return_value=native_scheme)
        native_scheme_cls.compute.return_value = 0
        result = Result()
        self.scheme.compute(inputs, native_scheme_cls, native_input, result)
        self.assertIsNotNone(self.scheme.inputs)
        self.assertIsNotNone(self.scheme.results)
        native_scheme_cls.assert_called_once_with(native_input)
        native_scheme.compute.assert_called_once()
        self.scheme.db_handler.insert_run.assert_called_once_with(inputs, result)

    @mock.patch("os.path.exists")
    @mock.patch("os.remove")
    @mock.patch("sys.exit")
    def test_check_status_and_clean_success(self, exit, remove, exists):
        exists.return_value = True
        self.scheme._check_status_and_clean(0, "recovery_file")
        exists.assert_called_once_with("recovery_file")
        remove.assert_called_once_with("recovery_file")
        exit.assert_not_called()

    @mock.patch("os.path.exists")
    @mock.patch("sys.exit")
    def test_check_status_and_clean_failure(self, exit, exists):
        exists.return_value = False
        self.scheme._check_status_and_clean(1, "recovery_file")
        exists.assert_not_called()
        exit.assert_called_once_with("Error while solving the dielectric theory")

    def test_save(self):
        self.scheme.results = Result()
        self.scheme._save()
        self.scheme.db_handler.insert_run.assert_called_once()

    def test_save_with_results(self):
        self.scheme.results = Result()
        self.scheme.inputs = Input(coupling=1.0, degeneracy=2.0)
        self.scheme._save()
        self.scheme.db_handler.insert_run.assert_called_once_with(
            self.scheme.inputs, self.scheme.results
        )

    def test_save_without_results(self):
        self.scheme.results = None
        self.scheme._save()
        self.scheme.db_handler.insert_run.assert_not_called()

    @mock.patch("qupled.base.Result.compute_rdf")
    def test_compute_rdf_with_default_grid(self, compute_rdf):
        self.scheme.db_handler = mock.MagicMock()
        self.scheme.results = Result()
        self.scheme.compute_rdf()
        compute_rdf.assert_called_once_with(None)
        self.scheme.db_handler.insert_results.assert_called_once_with(
            {"rdf": self.scheme.results.rdf, "rdf_grid": self.scheme.results.rdf_grid}
        )

    @mock.patch("qupled.base.Result.compute_rdf")
    def test_compute_rdf_with_custom_grid(self, compute_rdf):
        self.scheme.db_handler = mock.MagicMock()
        self.scheme.results = Result()
        rdf_grid = np.array([1, 2, 3])
        self.scheme.compute_rdf(rdf_grid)
        compute_rdf.assert_called_once_with(rdf_grid)
        self.scheme.db_handler.insert_results.assert_called_once_with(
            {"rdf": self.scheme.results.rdf, "rdf_grid": self.scheme.results.rdf_grid}
        )

    def test_compute_rdf_without_results(self):
        self.scheme.results = None
        self.scheme.compute_rdf()
        self.scheme.db_handler.insert_results.assert_not_called()


class TestIterativeScheme(TestCase):

    @mock.patch("qupled.util.DataBase.read_results")
    def test_get_initial_guess_with_default_database_name(self, read_results):
        run_id = "test_run"
        read_results.return_value = {
            "wvg": np.array([1, 2, 3]),
            "slfc": np.array([4, 5, 6]),
        }
        guess = IterativeScheme.get_initial_guess(run_id)
        self.assertTrue((guess.wvg == np.array([1, 2, 3])).all())
        self.assertTrue((guess.slfc == np.array([4, 5, 6])).all())
        read_results.assert_called_once_with(run_id, None, ["wvg", "slfc"])

    @mock.patch("qupled.util.DataBase.read_results")
    def test_get_initial_guess_with_custom_database_name(self, read_results):
        run_id = "test_run"
        database_name = "test_db"
        read_results.return_value = {
            "wvg": np.array([1, 2, 3]),
            "slfc": np.array([4, 5, 6]),
        }
        guess = IterativeScheme.get_initial_guess(run_id, database_name)
        self.assertTrue((guess.wvg == np.array([1, 2, 3])).all())
        self.assertTrue((guess.slfc == np.array([4, 5, 6])).all())
        read_results.assert_called_once_with(run_id, database_name, ["wvg", "slfc"])


class TestIterativeSchemeGuess(TestCase):

    def test_initialization(self):
        guess = IterativeScheme.Guess(wvg=np.array([1, 2, 3]), slfc=np.array([4, 5, 6]))
        self.assertTrue((guess.wvg == np.array([1, 2, 3])).all())
        self.assertTrue((guess.slfc == np.array([4, 5, 6])).all())

    def test_initialization_defaults(self):
        guess = IterativeScheme.Guess()
        self.assertIsNone(guess.wvg)
        self.assertIsNone(guess.slfc)

    @mock.patch("qupled.base.native.StlsGuess")
    def test_to_native(self, StlsGuess):
        native_guess = mock.MagicMock()
        StlsGuess.return_value = native_guess
        guess = IterativeScheme.Guess(wvg=np.array([1, 2, 3]), slfc=np.array([4, 5, 6]))
        result = guess.to_native()
        self.assertEqual(result, native_guess)
        self.assertTrue((result.wvg == np.array([1, 2, 3])).all())
        self.assertTrue((result.slfc == np.array([4, 5, 6])).all())


class TestQuantumIterativeSchemeGetInitialGuess(TestCase):

    @mock.patch("qupled.util.DataBase.read_run")
    def test_get_initial_guess_with_default_database_name(self, read_run):
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
        self.assertTrue((guess.wvg == np.array([1, 2, 3])).all())
        self.assertTrue((guess.ssf == np.array([4, 5, 6])).all())
        self.assertTrue((guess.adr == np.array([7, 8, 9])).all())
        self.assertEqual(guess.matsubara, 128)
        read_run.assert_called_once_with(
            run_id, None, ["matsubara"], ["wvg", "ssf", "adr"]
        )

    @mock.patch("qupled.util.DataBase.read_run")
    def test_get_initial_guess_with_custom_database_name(self, read_run):
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
        self.assertTrue((guess.wvg == np.array([1, 2, 3])).all())
        self.assertTrue((guess.ssf == np.array([4, 5, 6])).all())
        self.assertTrue((guess.adr == np.array([7, 8, 9])).all())
        self.assertEqual(guess.matsubara, 128)
        read_run.assert_called_once_with(
            run_id, database_name, ["matsubara"], ["wvg", "ssf", "adr"]
        )

class TestQuantumIterativeSchemeGuess(TestCase):

    def test_initialization(self):
        guess = QuantumIterativeScheme.Guess(
            wvg=np.array([1, 2, 3]),
            ssf=np.array([4, 5, 6]),
            adr=np.array([7, 8, 9]),
            matsubara=128,
        )
        self.assertTrue((guess.wvg == np.array([1, 2, 3])).all())
        self.assertTrue((guess.ssf == np.array([4, 5, 6])).all())
        self.assertTrue((guess.adr == np.array([7, 8, 9])).all())
        self.assertEqual(guess.matsubara, 128)

    def test_initialization_defaults(self):
        guess = QuantumIterativeScheme.Guess()
        self.assertIsNone(guess.wvg)
        self.assertIsNone(guess.ssf)
        self.assertIsNone(guess.adr)
        self.assertEqual(guess.matsubara, 0)

    @mock.patch("qupled.base.native.QstlsGuess")
    def test_to_native(self, QstlsGuess):
        native_guess = mock.MagicMock()
        QstlsGuess.return_value = native_guess
        guess = QuantumIterativeScheme.Guess(
            wvg=np.array([1, 2, 3]),
            ssf=np.array([4, 5, 6]),
            adr=np.array([7, 8, 9]),
            matsubara=128,
        )
        result = guess.to_native()
        self.assertEqual(result, native_guess)
        self.assertTrue((result.wvg == np.array([1, 2, 3])).all())
        self.assertTrue((result.ssf == np.array([4, 5, 6])).all())
        self.assertTrue((result.adr == np.array([7, 8, 9])).all())
        self.assertEqual(result.matsubara, 128)


if __name__ == "__main__":
    unittestmain()
