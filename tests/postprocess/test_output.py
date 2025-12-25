import pytest

from qupled.postprocess.output import DataBase


@pytest.fixture
def db_handler(mocker):
    yield mocker.patch("qupled.postprocess.output.DataBaseHandler")


@pytest.fixture
def scheme_tables(db_handler):
    return db_handler.return_value.scheme_tables


def test_inspect_runs(mocker, scheme_tables, db_handler):
    scheme_tables.inspect_runs.return_value = {"run1": "data1"}
    result = DataBase.inspect_runs(mocker.ANY, "test_db")
    assert result == {"run1": "data1"}
    db_handler.assert_called_once_with("test_db")


def test_read_run(mocker, scheme_tables, db_handler):
    scheme_tables.get_run.return_value = {
        "input1": "data1",
        "result1": "data2",
    }
    result = DataBase.read_run(mocker.ANY, 1, "test_db", ["input1"], ["result1"])
    assert result == {"input1": "data1", "result1": "data2"}
    db_handler.assert_called_once_with("test_db")


def test_read_inputs(mocker, scheme_tables, db_handler):
    scheme_tables.get_inputs.return_value = {"input1": "data1"}
    result = DataBase.read_inputs(mocker.ANY, 1, "test_db", ["input1"])
    assert result == {"input1": "data1"}
    db_handler.assert_called_once_with("test_db")


def test_read_results(mocker, scheme_tables, db_handler):
    scheme_tables.get_results.return_value = {"result1": "data1"}
    result = DataBase.read_results(mocker.ANY, 1, "test_db", ["result1"])
    assert result == {"result1": "data1"}
    db_handler.assert_called_once_with("test_db")
