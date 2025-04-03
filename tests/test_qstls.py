import numpy as np
from qupled.database import DataBaseHandler
import qupled.native as native
import qupled.qstls as qstls
import qupled.stls as stls


def test_qstls_inheritance():
    assert issubclass(qstls.Qstls, stls.Stls)


def test_qstls_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.__init__")
    scheme = qstls.Qstls()
    super_init.assert_called_once()
    assert isinstance(scheme.results, qstls.Result)
    assert scheme.native_scheme_cls == native.Qstls
    assert isinstance(scheme.native_inputs, native.QstlsInput)


def test_get_initial_guess_with_default_database_name(mocker):
    read_run = mocker.patch("qupled.output.DataBase.read_run")
    run_id = mocker.ANY
    read_run.return_value = {
        DataBaseHandler.INPUT_TABLE_NAME: {"matsubara": 128},
        DataBaseHandler.RESULT_TABLE_NAME: {
            "wvg": np.array([1, 2, 3]),
            "ssf": np.array([4, 5, 6]),
            "adr": np.array([7, 8, 9]),
        },
    }
    guess = qstls.Qstls.get_initial_guess(run_id)
    assert np.array_equal(guess.wvg, np.array([1, 2, 3]))
    assert np.array_equal(guess.ssf, np.array([4, 5, 6]))
    assert np.array_equal(guess.adr, np.array([7, 8, 9]))
    assert guess.matsubara == 128
    read_run.assert_called_once_with(run_id, None, ["matsubara"], ["wvg", "ssf", "adr"])


def test_get_initial_guess_with_custom_database_name(mocker):
    read_run = mocker.patch("qupled.output.DataBase.read_run")
    run_id = mocker.ANY
    database_name = mocker.ANY
    read_run.return_value = {
        DataBaseHandler.INPUT_TABLE_NAME: {"matsubara": 128},
        DataBaseHandler.RESULT_TABLE_NAME: {
            "wvg": np.array([1, 2, 3]),
            "ssf": np.array([4, 5, 6]),
            "adr": np.array([7, 8, 9]),
        },
    }
    guess = qstls.Qstls.get_initial_guess(run_id, database_name)
    assert np.array_equal(guess.wvg, np.array([1, 2, 3]))
    assert np.array_equal(guess.ssf, np.array([4, 5, 6]))
    assert np.array_equal(guess.adr, np.array([7, 8, 9]))
    assert guess.matsubara == 128
    read_run.assert_called_once_with(
        run_id, database_name, ["matsubara"], ["wvg", "ssf", "adr"]
    )


def test_qstls_input_inheritance():
    assert issubclass(qstls.Input, stls.Input)


def test_qstls_input_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    input = qstls.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "QSTLS"


def test_qstls_result_inheritance():
    assert issubclass(qstls.Result, stls.Result)


def test_qstls_result_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Result.__init__")
    result = qstls.Result()
    super_init.assert_called_once()
    assert result.adr is None
