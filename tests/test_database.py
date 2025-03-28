import os
import pytest
import sqlalchemy as sql
from sqlalchemy import inspect
from qupled.database import DataBaseHandler


@pytest.fixture
def db_handler():
    handler = DataBaseHandler()
    yield handler
    if os.path.exists(handler.database_name):
        os.remove(handler.database_name)


def test_database_handler_initialization_with_default_name(db_handler):
    assert db_handler.database_name == DataBaseHandler.DEFAULT_DATABASE_NAME
    assert db_handler.engine.url.database == DataBaseHandler.DEFAULT_DATABASE_NAME
    assert db_handler.run_id is None
    inspector = inspect(db_handler.engine)
    assert set(inspector.get_table_names()) == {
        DataBaseHandler.RUNS_TABLE_NAME,
        DataBaseHandler.INPUTS_TABLE_NAME,
        DataBaseHandler.RESULTS_TABLE_NAME,
    }


def test_database_handler_initialization_with_custom_name():
    database_name = "custom.db"
    db_handler = DataBaseHandler(database_name="custom.db")
    assert db_handler.database_name == database_name
    assert db_handler.engine.url.database == database_name
    assert db_handler.run_id is None
    inspector = inspect(db_handler.engine)
    assert set(inspector.get_table_names()) == {
        DataBaseHandler.RUNS_TABLE_NAME,
        DataBaseHandler.INPUTS_TABLE_NAME,
        DataBaseHandler.RESULTS_TABLE_NAME,
    }
    if os.path.exists(db_handler.database_name):
        os.remove(db_handler.database_name)


def test_set_sqlite_pragma_valid_engine():
    engine = sql.create_engine("sqlite:///:memory:")
    DataBaseHandler._set_sqlite_pragma(engine)
    with engine.connect() as connection:
        result = connection.execute(sql.text("PRAGMA foreign_keys")).fetchone()
        assert result[0] == 1


def test_set_sqlite_pragma_invalid_engine():
    invalid_engine = None
    with pytest.raises(sql.exc.InvalidRequestError):
        DataBaseHandler._set_sqlite_pragma(invalid_engine)


def test_insert_run(mocker, db_handler):
    insert_run = mocker.patch.object(db_handler, "_insert_run")
    insert_inputs = mocker.patch.object(db_handler, "insert_inputs")
    insert_results = mocker.patch.object(db_handler, "insert_results")
    inputs = mocker.MagicMock()
    results = mocker.MagicMock()
    db_handler.insert_run(inputs, results)
    insert_run.assert_called_once_with(inputs)
    insert_inputs.assert_called_once_with(inputs.__dict__)
    insert_results.assert_called_once_with(results.__dict__)


def test_insert_inputs(mocker, db_handler):
    db_handler.input_table = mocker.MagicMock()
    insert_from_dict = mocker.patch.object(db_handler, "_insert_from_dict")
    mocker.patch.object(db_handler, "_to_json", side_effect=lambda x: f"json({x})")
    inputs = mocker.MagicMock()
    db_handler.insert_inputs(inputs)
    insert_from_dict.assert_called_once()
    called_table, called_inputs, called_mapper = insert_from_dict.call_args[0]
    assert called_table == db_handler.result_table
    assert called_inputs == inputs
    assert callable(called_mapper)
    assert called_mapper("value1") == "json(value1)"
    assert called_mapper(123) == "json(123)"


def test_insert_results(mocker, db_handler):
    db_handler.result_table = mocker.MagicMock()
    insert_from_dict = mocker.patch.object(db_handler, "_insert_from_dict")
    mocker.patch.object(db_handler, "_to_bytes", side_effect=lambda x: f"bytes({x})")
    results = mocker.MagicMock()
    db_handler.insert_results(results)
    insert_from_dict.assert_called_once()
    called_table, called_inputs, called_mapper = insert_from_dict.call_args[0]
    assert called_table == db_handler.result_table
    assert called_inputs == results
    assert callable(called_mapper)
    assert called_mapper("value1") == "bytes(value1)"
    assert called_mapper(123) == "bytes(123)"


def test_inspect_runs(mocker, db_handler):
    sql_select = mocker.patch("sqlalchemy.select")
    mock_runs = [{"id": 1, "status": "done"}, {"id": 2, "status": "failed"}]
    db_handler.run_table = mocker.MagicMock()
    execute = mocker.patch.object(db_handler, "_execute")
    mock_result = mocker.MagicMock()
    mock_result.mappings.return_value.all.return_value = mock_runs
    execute.return_value = mock_result
    runs = db_handler.inspect_runs()
    sql_select.assert_called_once_with(db_handler.run_table)
    mock_result.mappings.return_value.all.assert_called_once()
    assert runs == mock_runs


def test_get_run_with_existing_run(mocker, db_handler):
    run_id = 1
    sql_select = mocker.patch("sqlalchemy.select")
    execute = mocker.patch.object(db_handler, "_execute")
    db_handler.run_table = mocker.MagicMock()
    mock_statement = mocker.MagicMock()
    sql_select.return_value.where.return_value = mock_statement
    mock_result = mocker.MagicMock()
    mock_result.mappings.return_value.first.return_value = {"key": "value"}
    execute.return_value = mock_result
    inputs = mocker.MagicMock()
    get_inputs = mocker.patch.object(db_handler, "get_inputs", return_value=inputs)
    results = mocker.MagicMock()
    get_results = mocker.patch.object(db_handler, "get_results", return_value=results)
    run = db_handler.get_run(run_id, None, None)
    sql_select.assert_called_once_with(db_handler.run_table)
    get_inputs.assert_called_once_with(run_id, names=None)
    get_results.assert_called_once_with(run_id, names=None)
    execute.assert_called_once_with(mock_statement)
    mock_result.mappings.return_value.first.assert_called_once()
    assert run == {
        DataBaseHandler.RUNS_TABLE_NAME: {"key": "value"},
        DataBaseHandler.INPUTS_TABLE_NAME: inputs,
        DataBaseHandler.RESULTS_TABLE_NAME: results,
    }


def test_get_run_with_non_existing_run(mocker, db_handler):
    run_id = 1
    sql_select = mocker.patch("sqlalchemy.select")
    db_handler.run_table = mocker.MagicMock()
    mock_statement = mocker.MagicMock()
    sql_select.return_value.where.return_value = mock_statement
    execute = mocker.patch.object(db_handler, "_execute")
    mock_result = mocker.MagicMock()
    mock_result.mappings.return_value.first.return_value = None
    execute.return_value = mock_result
    run = db_handler.get_run(run_id, None, None)
    sql_select.assert_called_once_with(db_handler.run_table)
    execute.assert_called_once_with(mock_statement)
    mock_result.mappings.return_value.first.assert_called_once()
    assert run == {}


def test_get_inputs(mocker, db_handler):
    run_id = 1
    names = ["name"]
    db_handler.input_table = mocker.MagicMock()
    expected_inputs = {"key": "value"}
    get = mocker.patch.object(db_handler, "_get", return_value=expected_inputs)
    mocker.patch.object(
        db_handler, "_from_json", side_effect=lambda x: f"from_json({x})"
    )
    inputs = db_handler.get_inputs(run_id, names)
    get.assert_called_once()
    called_table, called_run_id, called_names, called_mapper = get.call_args[0]
    assert called_table == db_handler.input_table
    assert called_run_id == run_id
    assert called_names == names
    assert callable(called_mapper)
    assert called_mapper("value1") == "from_json(value1)"
    assert called_mapper(123) == "from_json(123)"
    assert inputs == expected_inputs


def test_get_results(mocker, db_handler):
    run_id = 1
    names = ["name"]
    db_handler.result_table = mocker.MagicMock()
    expected_results = {"key": "value"}
    get = mocker.patch.object(db_handler, "_get", return_value=expected_results)
    mocker.patch.object(
        db_handler, "_from_bytes", side_effect=lambda x: f"from_bytes({x})"
    )
    results = db_handler.get_results(run_id, names)
    get.assert_called_once()
    called_table, called_run_id, called_names, called_mapper = get.call_args[0]
    assert called_table == db_handler.result_table
    assert called_run_id == run_id
    assert called_names == names
    assert callable(called_mapper)
    assert called_mapper("value1") == "from_bytes(value1)"
    assert called_mapper(123) == "from_bytes(123)"
    assert results == expected_results


def test_delete_run_with_existing_run(mocker, db_handler):
    run_id = 1
    db_handler.run_table = mocker.MagicMock()
    db_handler.run_table.result_value.c = mocker.MagicMock()
    statement = mocker.MagicMock()
    sql_delete = mocker.patch("sqlalchemy.delete")
    sql_delete.return_value.where.return_value = statement
    execute = mocker.patch.object(db_handler, "_execute")
    db_handler.delete_run(run_id)
    sql_delete.assert_called_once_with(db_handler.run_table)
    execute.assert_called_once_with(statement)
