import pytest
from qupled.output import DataBase


@pytest.fixture
def db_handler(mocker):
    yield mocker.patch("qupled.database.DataBaseHandler")


def test_inspect_runs(db_handler):
    db_handler.return_value.inspect_runs.return_value = {"run1": "data1"}
    result = DataBase.inspect_runs("test_db")
    assert result == {"run1": "data1"}
    db_handler.assert_called_once_with("test_db")


def test_read_run(db_handler):
    db_handler.return_value.get_run.return_value = {
        "input1": "data1",
        "result1": "data2",
    }
    result = DataBase.read_run(1, "test_db", ["input1"], ["result1"])
    assert result == {"input1": "data1", "result1": "data2"}
    db_handler.assert_called_once_with("test_db")


def test_read_inputs(db_handler):
    db_handler.return_value.get_inputs.return_value = {"input1": "data1"}
    result = DataBase.read_inputs(1, "test_db", ["input1"])
    assert result == {"input1": "data1"}
    db_handler.assert_called_once_with("test_db")


def test_read_results(db_handler):
    db_handler.return_value.get_results.return_value = {"result1": "data1"}
    result = DataBase.read_results(1, "test_db", ["result1"])
    assert result == {"result1": "data1"}
    db_handler.assert_called_once_with("test_db")
