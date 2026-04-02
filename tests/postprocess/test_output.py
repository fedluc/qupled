import pytest

from qupled.postprocess.output import DataBase, OutputType


@pytest.fixture
def db_handler(mocker):
    yield mocker.patch("qupled.postprocess.output.DataBaseHandler")


@pytest.fixture
def scheme_tables(db_handler):
    return db_handler.return_value.scheme_tables


@pytest.fixture
def fsc_tables(db_handler):
    return db_handler.return_value.fsc_tables


@pytest.mark.unit
def test_inspect_runs_scheme(scheme_tables, db_handler):
    scheme_tables.inspect_runs.return_value = {"run1": "data1"}
    result = DataBase.inspect_runs(OutputType.SCHEME, "test_db")
    assert result == {"run1": "data1"}
    db_handler.assert_called_once_with("test_db")
    scheme_tables.inspect_runs.assert_called_once()


@pytest.mark.unit
def test_inspect_runs_fsc(fsc_tables, db_handler):
    fsc_tables.inspect_runs.return_value = {"run1": "data1"}
    result = DataBase.inspect_runs(OutputType.FINITE_SIZE_CORRECTION, "test_db")
    assert result == {"run1": "data1"}
    db_handler.assert_called_once_with("test_db")
    fsc_tables.inspect_runs.assert_called_once()


@pytest.mark.unit
def test_read_run_scheme(scheme_tables, db_handler):
    scheme_tables.get_run.return_value = {"input1": "data1", "result1": "data2"}
    result = DataBase.read_run(1, OutputType.SCHEME, "test_db", ["input1"], ["result1"])
    assert result == {"input1": "data1", "result1": "data2"}
    db_handler.assert_called_once_with("test_db")
    scheme_tables.get_run.assert_called_once_with(1, ["input1"], ["result1"])


@pytest.mark.unit
def test_read_run_fsc(fsc_tables, db_handler):
    fsc_tables.get_run.return_value = {"input1": "data1", "result1": "data2"}
    result = DataBase.read_run(
        1, OutputType.FINITE_SIZE_CORRECTION, "test_db", ["input1"], ["result1"]
    )
    assert result == {"input1": "data1", "result1": "data2"}
    db_handler.assert_called_once_with("test_db")
    fsc_tables.get_run.assert_called_once_with(1, ["input1"], ["result1"])


@pytest.mark.unit
def test_read_inputs_scheme(scheme_tables, db_handler):
    scheme_tables.get_inputs.return_value = {"input1": "data1"}
    result = DataBase.read_inputs(1, OutputType.SCHEME, "test_db", ["input1"])
    assert result == {"input1": "data1"}
    db_handler.assert_called_once_with("test_db")
    scheme_tables.get_inputs.assert_called_once_with(1, ["input1"])


@pytest.mark.unit
def test_read_inputs_fsc(fsc_tables, db_handler):
    fsc_tables.get_inputs.return_value = {"input1": "data1"}
    result = DataBase.read_inputs(
        1, OutputType.FINITE_SIZE_CORRECTION, "test_db", ["input1"]
    )
    assert result == {"input1": "data1"}
    db_handler.assert_called_once_with("test_db")
    fsc_tables.get_inputs.assert_called_once_with(1, ["input1"])


@pytest.mark.unit
def test_read_results_scheme(scheme_tables, db_handler):
    scheme_tables.get_results.return_value = {"result1": "data1"}
    result = DataBase.read_results(1, OutputType.SCHEME, "test_db", ["result1"])
    assert result == {"result1": "data1"}
    db_handler.assert_called_once_with("test_db")
    scheme_tables.get_results.assert_called_once_with(1, ["result1"])


@pytest.mark.unit
def test_read_results_fsc(fsc_tables, db_handler):
    fsc_tables.get_results.return_value = {"result1": "data1"}
    result = DataBase.read_results(
        1, OutputType.FINITE_SIZE_CORRECTION, "test_db", ["result1"]
    )
    assert result == {"result1": "data1"}
    db_handler.assert_called_once_with("test_db")
    fsc_tables.get_results.assert_called_once_with(1, ["result1"])


@pytest.mark.unit
def test_delete_run_scheme(scheme_tables, db_handler):
    DataBase.delete_run(1, OutputType.SCHEME, "test_db")
    db_handler.assert_called_once_with("test_db")
    scheme_tables.delete_run.assert_called_once_with(1)


@pytest.mark.unit
def test_delete_run_fsc(fsc_tables, db_handler):
    DataBase.delete_run(1, OutputType.FINITE_SIZE_CORRECTION, "test_db")
    db_handler.assert_called_once_with("test_db")
    fsc_tables.delete_run.assert_called_once_with(1)
