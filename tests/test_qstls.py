import pytest

import qupled.native as native
import qupled.qstls as qstls
import qupled.stls as stls
from qupled.base_tables import TableKeys


@pytest.fixture
def scheme():
    scheme = qstls.Solver()
    return scheme


def test_qstls_inheritance():
    assert issubclass(qstls.Solver, stls.Solver)


def test_qstls_initialization(mocker, scheme):
    super_init = mocker.patch("qupled.stls.Solver.__init__")
    scheme = qstls.Solver()
    super_init.assert_called_once()
    assert isinstance(scheme.results, stls.Result)
    assert scheme.native_scheme_cls == native.Qstls
    assert scheme.native_inputs_cls == native.QstlsInput


def test_compute(mocker, scheme):
    find_fixed_adr_in_database = mocker.patch(
        "qupled.qstls.Solver.find_fixed_adr_in_database"
    )
    super_compute = mocker.patch("qupled.stls.Solver.compute")
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
    db_handler_mock.inspect_scheme_runs.return_value = [
        {
            TableKeys.DEGENERACY.value: 2.0,
            TableKeys.THEORY.value: mocker.ANY,
            TableKeys.PRIMARY_KEY.value: 1,
        }
    ]
    db_handler_mock.get_scheme_inputs.return_value = {
        "cutoff": 10,
        "matsubara": 128,
        "resolution": 0.01,
    }
    scheme.find_fixed_adr_in_database(inputs)
    assert inputs.fixed_run_id == 1
    db_handler_mock.inspect_scheme_runs.assert_called_once()
    db_handler_mock.get_scheme_inputs.assert_called_once_with(1)


def test_find_fixed_adr_in_database_no_match(mocker, scheme):
    db_handler_mock = mocker.Mock()
    scheme.db_handler = db_handler_mock
    inputs = qstls.Input(coupling=mocker.ANY, degeneracy=2.0)
    inputs.cutoff = 10
    inputs.matsubara = 128
    inputs.resolution = 0.01
    db_handler_mock.inspect_scheme_runs.return_value = [
        {
            TableKeys.DEGENERACY.value: 3.0,
            TableKeys.THEORY.value: mocker.ANY,
            TableKeys.PRIMARY_KEY.value: mocker.ANY,
        }
    ]
    scheme.find_fixed_adr_in_database(inputs)
    assert inputs.fixed_run_id is None
    db_handler_mock.inspect_scheme_runs.assert_called_once()
    db_handler_mock.get_scheme_inputs.assert_not_called()


def test_qstls_input_inheritance():
    assert issubclass(qstls.Input, stls.Input)


def test_qstls_input_initialization(mocker):
    input = qstls.Input(mocker.ANY, mocker.ANY)
    assert input.theory == "QSTLS"
