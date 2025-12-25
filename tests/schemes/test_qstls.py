import pytest

from qupled import native
from qupled.database.scheme_tables import TableKeys, BaseTableKeys
from qupled.schemes import qstls, stls


@pytest.fixture
def scheme():
    scheme = qstls.Solver()
    return scheme


def test_qstls_inheritance():
    assert issubclass(qstls.Solver, stls.Solver)


def test_qstls_initialization(mocker, scheme):
    super_init = mocker.patch("qupled.schemes.stls.Solver.__init__")
    scheme = qstls.Solver()
    super_init.assert_called_once()
    assert isinstance(scheme.results, stls.Result)
    assert scheme.native_scheme_cls == native.Qstls
    assert scheme.native_inputs_cls == native.QstlsInput


def test_compute(mocker, scheme):
    find_fixed_adr_in_database = mocker.patch(
        "qupled.schemes.qstls.Solver.find_fixed_adr_in_database"
    )
    super_compute = mocker.patch("qupled.schemes.stls.Solver.compute")
    inputs = mocker.ANY
    scheme.compute(inputs)
    find_fixed_adr_in_database.assert_called_once_with(inputs)
    super_compute.assert_called_once_with(inputs)


def test_find_fixed_adr_in_database_match_found(mocker, scheme):
    db_handler_mock = mocker.Mock()
    scheme_tables = db_handler_mock.scheme_tables
    scheme.db_handler = db_handler_mock
    inputs = qstls.Input(coupling=mocker.ANY, degeneracy=2.0)
    inputs.cutoff = 10
    inputs.matsubara = 128
    inputs.resolution = 0.01
    scheme_tables.inspect_runs.return_value = [
        {
            TableKeys.DEGENERACY.value: 2.0,
            TableKeys.THEORY.value: mocker.ANY,
            BaseTableKeys.PRIMARY_KEY.value: 1,
        }
    ]
    scheme_tables.get_inputs.return_value = {
        "cutoff": inputs.cutoff,
        "matsubara": inputs.matsubara,
        "resolution": inputs.resolution,
    }
    scheme.find_fixed_adr_in_database(inputs)
    assert inputs.fixed_run_id == 1
    scheme_tables.inspect_runs.assert_called_once()
    scheme_tables.get_inputs.assert_called_once_with(1)


def test_find_fixed_adr_in_database_no_match(mocker, scheme):
    db_handler_mock = mocker.Mock()
    scheme_tables = db_handler_mock.scheme_tables
    scheme.db_handler = db_handler_mock
    inputs = mocker.Mock()
    scheme_tables.inspect_runs.return_value = [
        {
            TableKeys.DEGENERACY.value: 3.0,
            TableKeys.THEORY.value: mocker.ANY,
            BaseTableKeys.PRIMARY_KEY.value: mocker.ANY,
        }
    ]
    scheme.find_fixed_adr_in_database(inputs)
    assert inputs.fixed_run_id is None
    scheme_tables.inspect_runs.assert_called_once()
    scheme_tables.get_inputs.assert_not_called()


def test_qstls_input_inheritance():
    assert issubclass(qstls.Input, stls.Input)


def test_qstls_input_initialization(mocker):
    input = qstls.Input(mocker.ANY, mocker.ANY)
    assert input.theory == "QSTLS"
