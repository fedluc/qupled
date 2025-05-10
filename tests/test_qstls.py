import pytest

import qupled.native as native
import qupled.qstls as qstls
import qupled.stls as stls
from qupled.database import DataBaseHandler


@pytest.fixture
def scheme():
    scheme = qstls.Qstls()
    return scheme


def test_qstls_inheritance():
    assert issubclass(qstls.Qstls, stls.Stls)


def test_qstls_initialization(mocker, scheme):
    super_init = mocker.patch("qupled.stls.Stls.__init__")
    scheme = qstls.Qstls()
    super_init.assert_called_once()
    assert isinstance(scheme.results, stls.Result)
    assert scheme.native_scheme_cls == native.Qstls
    assert scheme.native_inputs_cls == native.QstlsInput


def test_compute(mocker, scheme):
    find_fixed_adr_in_database = mocker.patch(
        "qupled.qstls.Qstls.find_fixed_adr_in_database"
    )
    super_compute = mocker.patch("qupled.stls.Stls.compute")
    inputs = mocker.ANY
    scheme.compute(inputs)
    find_fixed_adr_in_database.assert_called_once_with(inputs)
    super_compute.assert_called_once_with(inputs)


def test_find_fixed_adr_in_database_match_found(mocker, scheme):
    db_handler_mock = mocker.Mock()
    scheme.db_handler = db_handler_mock
    inputs = qstls.Input(coupling=mocker.ANY, degeneracy=2.0)
    inputs.cutoff = 10
    inputs.matsubara = 128
    inputs.resolution = 0.01
    database_keys = DataBaseHandler.TableKeys
    db_handler_mock.inspect_runs.return_value = [
        {
            database_keys.DEGENERACY.value: 2.0,
            database_keys.THEORY.value: mocker.ANY,
            database_keys.PRIMARY_KEY.value: 1,
        }
    ]
    db_handler_mock.get_inputs.return_value = {
        "cutoff": 10,
        "matsubara": 128,
        "resolution": 0.01,
    }
    scheme.find_fixed_adr_in_database(inputs)
    assert inputs.fixed_run_id == 1
    db_handler_mock.inspect_runs.assert_called_once()
    db_handler_mock.get_inputs.assert_called_once_with(1)


def test_find_fixed_adr_in_database_no_match(mocker, scheme):
    db_handler_mock = mocker.Mock()
    scheme.db_handler = db_handler_mock
    inputs = qstls.Input(coupling=mocker.ANY, degeneracy=2.0)
    inputs.cutoff = 10
    inputs.matsubara = 128
    inputs.resolution = 0.01
    database_keys = DataBaseHandler.TableKeys
    db_handler_mock.inspect_runs.return_value = [
        {
            database_keys.DEGENERACY.value: 3.0,
            database_keys.THEORY.value: mocker.ANY,
            database_keys.PRIMARY_KEY.value: mocker.ANY,
        }
    ]
    scheme.find_fixed_adr_in_database(inputs)
    assert inputs.fixed_run_id is None
    db_handler_mock.inspect_runs.assert_called_once()
    db_handler_mock.get_inputs.assert_not_called()


def test_qstls_input_inheritance():
    assert issubclass(qstls.Input, stls.Input)


def test_qstls_input_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    input = qstls.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "QSTLS"
