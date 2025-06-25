import datetime
import io
import json
import os
import struct

import blosc2
import numpy as np
import pytest
import sqlalchemy as sql
from sqlalchemy import inspect

from qupled.database import DataBaseHandler

# Unit tests


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
        DataBaseHandler.RUN_TABLE_NAME,
        DataBaseHandler.INPUT_TABLE_NAME,
        DataBaseHandler.RESULT_TABLE_NAME,
    }


def test_database_handler_initialization_with_custom_name():
    database_name = "custom.db"
    db_handler = DataBaseHandler(database_name="custom.db")
    assert db_handler.database_name == database_name
    assert db_handler.engine.url.database == database_name
    assert db_handler.run_id is None
    inspector = inspect(db_handler.engine)
    assert set(inspector.get_table_names()) == {
        DataBaseHandler.RUN_TABLE_NAME,
        DataBaseHandler.INPUT_TABLE_NAME,
        DataBaseHandler.RESULT_TABLE_NAME,
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
    inputs = mocker.ANY
    db_handler.insert_run(inputs)
    insert_run.assert_called_once_with(inputs, DataBaseHandler.RunStatus.RUNNING)
    insert_inputs.assert_called_once_with(inputs.__dict__)


def test_insert_inputs(mocker, db_handler):
    db_handler.run_id = mocker.ANY
    insert_from_dict = mocker.patch.object(db_handler, "_insert_from_dict")
    mocker.patch.object(db_handler, "_to_json", side_effect=lambda x: f"json({x})")
    inputs = mocker.ANY
    db_handler.insert_inputs(inputs)
    insert_from_dict.assert_called_once()
    called_table, called_inputs, called_mapper = insert_from_dict.call_args[0]
    assert called_table == db_handler.input_table
    assert called_inputs == inputs
    assert callable(called_mapper)
    assert called_mapper("value1") == "json(value1)"
    assert called_mapper(123) == "json(123)"


def test_insert_inputs_without_run_id(mocker, db_handler):
    insert_from_dict = mocker.patch.object(db_handler, "_insert_from_dict")
    mocker.patch.object(db_handler, "_to_json", side_effect=lambda x: f"json({x})")
    db_handler.insert_inputs(mocker.ANY)
    insert_from_dict.assert_not_called()


def test_insert_results(mocker, db_handler):
    db_handler.run_id = mocker.ANY
    insert_from_dict = mocker.patch.object(db_handler, "_insert_from_dict")
    mocker.patch.object(db_handler, "_to_bytes", side_effect=lambda x: f"bytes({x})")
    results = mocker.ANY
    db_handler.insert_results(results)
    insert_from_dict.assert_called_once()
    called_table, called_results, called_mapper, called_conflict_mode = (
        insert_from_dict.call_args[0]
    )
    assert called_table == db_handler.result_table
    assert called_results == results
    assert callable(called_mapper)
    assert called_mapper("value1") == "bytes(value1)"
    assert called_mapper(123) == "bytes(123)"
    assert called_conflict_mode == DataBaseHandler.ConflictMode.FAIL


def test_insert_results_without_run_id(mocker, db_handler):
    insert_from_dict = mocker.patch.object(db_handler, "_insert_from_dict")
    mocker.patch.object(db_handler, "_to_json", side_effect=lambda x: f"json({x})")
    db_handler.insert_results(mocker.ANY)
    insert_from_dict.assert_not_called()


def test_inspect_runs(mocker, db_handler):
    sql_select = mocker.patch("sqlalchemy.select")
    mock_runs = [{"id": 1, "status": "done"}, {"id": 2, "status": "failed"}]
    execute = mocker.patch.object(db_handler, "_execute")
    mock_result = mocker.Mock()
    mock_result.mappings.return_value.all.return_value = mock_runs
    execute.return_value = mock_result
    runs = db_handler.inspect_runs()
    sql_select.assert_called_once_with(db_handler.run_table)
    mock_result.mappings.return_value.all.assert_called_once()
    assert runs == mock_runs


def test_update_run_status_with_run_id(mocker, db_handler):
    db_handler.run_id = 1
    status = DataBaseHandler.RunStatus.SUCCESS
    sql_update = mocker.patch("sqlalchemy.update")
    statement = sql_update.return_value.where.return_value.values.return_value
    execute = mocker.patch.object(db_handler, "_execute")
    db_handler.update_run_status(status)
    sql_update.assert_called_once_with(db_handler.run_table)
    sql_update.return_value.where.return_value.values.assert_called_once_with(
        {DataBaseHandler.TableKeys.STATUS.value: status.value}
    )
    execute.assert_called_once_with(statement)


def test_update_run_status_without_run_id(mocker, db_handler):
    db_handler.run_id = None
    status = 0
    sql_update = mocker.patch("sqlalchemy.update")
    execute = mocker.patch.object(db_handler, "_execute")
    db_handler.update_run_status(status)
    sql_update.assert_not_called()
    execute.assert_not_called()


def test_get_run_with_existing_run(mocker, db_handler):
    run_id = 1
    sql_select = mocker.patch("sqlalchemy.select")
    execute = mocker.patch.object(db_handler, "_execute")
    statement = sql_select.return_value.where.return_value
    mock_result = mocker.Mock()
    mock_result.mappings.return_value.first.return_value = {"key": "value"}
    execute.return_value = mock_result
    inputs = mocker.ANY
    get_inputs = mocker.patch.object(db_handler, "get_inputs", return_value=inputs)
    results = mocker.ANY
    get_results = mocker.patch.object(db_handler, "get_results", return_value=results)
    run = db_handler.get_run(run_id, None, None)
    sql_select.assert_called_once_with(db_handler.run_table)
    get_inputs.assert_called_once_with(run_id, names=None)
    get_results.assert_called_once_with(run_id, names=None)
    execute.assert_called_once_with(statement)
    mock_result.mappings.return_value.first.assert_called_once()
    assert run == {
        DataBaseHandler.RUN_TABLE_NAME: {"key": "value"},
        DataBaseHandler.INPUT_TABLE_NAME: inputs,
        DataBaseHandler.RESULT_TABLE_NAME: results,
    }


def test_get_run_with_non_existing_run(mocker, db_handler):
    run_id = 1
    sql_select = mocker.patch("sqlalchemy.select")
    statement = sql_select.return_value.where.return_value
    execute = mocker.patch.object(db_handler, "_execute")
    mock_result = mocker.Mock()
    mock_result.mappings.return_value.first.return_value = None
    execute.return_value = mock_result
    run = db_handler.get_run(run_id, None, None)
    sql_select.assert_called_once_with(db_handler.run_table)
    execute.assert_called_once_with(statement)
    mock_result.mappings.return_value.first.assert_called_once()
    assert run == {}


def test_get_inputs(mocker, db_handler):
    run_id = 1
    names = ["name"]
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


def test_delete_run(mocker, db_handler):
    run_id = 1
    db_handler.run_table = mocker.MagicMock()
    db_handler.run_table.result_value.c = mocker.ANY
    sql_delete = mocker.patch("sqlalchemy.delete")
    statement = sql_delete.return_value.where.return_value
    execute = mocker.patch.object(db_handler, "_execute")
    db_handler.delete_run(run_id)
    sql_delete.assert_called_once_with(db_handler.run_table)
    execute.assert_called_once_with(statement)


def test_build_run_table_columns(db_handler):
    table = db_handler.run_table
    columns = {col.name for col in table.columns}
    expected_columns = {
        DataBaseHandler.TableKeys.PRIMARY_KEY.value,
        DataBaseHandler.TableKeys.THEORY.value,
        DataBaseHandler.TableKeys.COUPLING.value,
        DataBaseHandler.TableKeys.DEGENERACY.value,
        DataBaseHandler.TableKeys.DATE.value,
        DataBaseHandler.TableKeys.TIME.value,
        DataBaseHandler.TableKeys.STATUS.value,
    }
    assert columns == expected_columns


def test_build_run_table_primary_key(db_handler):
    table = db_handler.run_table
    primary_key_columns = {col.name for col in table.primary_key.columns}
    assert primary_key_columns == {DataBaseHandler.TableKeys.PRIMARY_KEY.value}


def test_build_run_table_column_types(db_handler):
    table = db_handler.run_table
    assert isinstance(
        table.c[DataBaseHandler.TableKeys.PRIMARY_KEY.value].type, sql.Integer
    )
    assert isinstance(table.c[DataBaseHandler.TableKeys.THEORY.value].type, sql.String)
    assert isinstance(table.c[DataBaseHandler.TableKeys.COUPLING.value].type, sql.Float)
    assert isinstance(
        table.c[DataBaseHandler.TableKeys.DEGENERACY.value].type, sql.Float
    )
    assert isinstance(table.c[DataBaseHandler.TableKeys.DATE.value].type, sql.String)
    assert isinstance(table.c[DataBaseHandler.TableKeys.TIME.value].type, sql.String)
    assert isinstance(table.c[DataBaseHandler.TableKeys.STATUS.value].type, sql.String)


def test_build_inputs_table_columns(db_handler):
    table = db_handler.input_table
    columns = {col.name for col in table.columns}
    expected_columns = {
        DataBaseHandler.TableKeys.RUN_ID.value,
        DataBaseHandler.TableKeys.NAME.value,
        DataBaseHandler.TableKeys.VALUE.value,
    }
    assert columns == expected_columns


def test_build_inputs_table_primary_key(db_handler):
    table = db_handler.input_table
    primary_key_columns = {col.name for col in table.primary_key.columns}
    assert primary_key_columns == {
        DataBaseHandler.TableKeys.RUN_ID.value,
        DataBaseHandler.TableKeys.NAME.value,
    }


def test_build_inputs_table_column_types(db_handler):
    table = db_handler.input_table
    assert isinstance(table.c[DataBaseHandler.TableKeys.RUN_ID.value].type, sql.Integer)
    assert isinstance(table.c[DataBaseHandler.TableKeys.NAME.value].type, sql.String)
    assert isinstance(table.c[DataBaseHandler.TableKeys.VALUE.value].type, sql.JSON)


def test_build_results_table_columns(db_handler):
    table = db_handler.result_table
    columns = {col.name for col in table.columns}
    expected_columns = {
        DataBaseHandler.TableKeys.RUN_ID.value,
        DataBaseHandler.TableKeys.NAME.value,
        DataBaseHandler.TableKeys.VALUE.value,
    }
    assert columns == expected_columns


def test_build_results_table_primary_key(db_handler):
    table = db_handler.result_table
    primary_key_columns = {col.name for col in table.primary_key.columns}
    assert primary_key_columns == {
        DataBaseHandler.TableKeys.RUN_ID.value,
        DataBaseHandler.TableKeys.NAME.value,
    }


def test_build_results_table_column_types(db_handler):
    table = db_handler.result_table
    assert isinstance(table.c[DataBaseHandler.TableKeys.RUN_ID.value].type, sql.Integer)
    assert isinstance(table.c[DataBaseHandler.TableKeys.NAME.value].type, sql.String)
    assert isinstance(
        table.c[DataBaseHandler.TableKeys.VALUE.value].type, sql.LargeBinary
    )


def test_insert_run(mocker, db_handler):
    run_id = 42
    fixed_now = datetime.datetime(2023, 1, 1, 12, 0, 0)
    mock_datetime = mocker.patch("qupled.database.datetime")
    mock_datetime.now.return_value = fixed_now
    mock_datetime.side_effect = lambda *args, **kwargs: datetime.datetime(
        *args, **kwargs
    )
    inputs = mocker.Mock()
    status = mocker.Mock()
    data = mocker.ANY
    sql_insert = mocker.patch("sqlalchemy.insert")
    statement = sql_insert.return_value.values.return_value
    result = mocker.Mock()
    execute = mocker.patch.object(db_handler, "_execute", return_value=result)
    result.inserted_primary_key = [run_id]
    db_handler._insert_run(inputs, status)
    mock_datetime.now.assert_called_once()
    sql_insert.assert_called_once_with(db_handler.run_table)
    sql_insert.return_value.values.assert_called_once_with(data)
    execute.assert_called_once_with(statement)
    assert db_handler.run_id == run_id


def test_insert_from_dict_with_valid_data(mocker, db_handler):
    table = mocker.ANY
    data = {"key1": "value1", "key2": "value2"}
    sql_mapping = mocker.Mock(side_effect=lambda x: f"mapped({x})")
    insert = mocker.patch.object(db_handler, "_insert")
    db_handler._insert_from_dict(table, data, sql_mapping)
    sql_mapping.assert_has_calls([mocker.call("value1"), mocker.call("value2")])
    insert.assert_has_calls(
        [
            mocker.call(table, "key1", "mapped(value1)", mocker.ANY),
            mocker.call(table, "key2", "mapped(value2)", mocker.ANY),
        ]
    )
    assert insert.call_count == 2


def test_insert_from_dict_with_empty_data(mocker, db_handler):
    data = {}
    sql_mapping = mocker.Mock()
    insert = mocker.patch.object(db_handler, "_insert")
    db_handler._insert_from_dict(mocker.ANY, data, sql_mapping)
    sql_mapping.assert_not_called()
    insert.assert_not_called()


def test_insert_with_conflict_mode_fail(mocker, db_handler):
    table = mocker.ANY
    name = "test_name"
    value = "test_value"
    run_id = 1
    db_handler.run_id = run_id
    sqlite_insert = mocker.patch("qupled.database.sqlite_insert")
    statement = sqlite_insert.return_value.values.return_value
    execute = mocker.patch.object(db_handler, "_execute")
    db_handler._insert(table, name, value)
    sqlite_insert.assert_called_once_with(table)
    sqlite_insert.return_value.values.assert_called_once_with(
        {
            db_handler.TableKeys.RUN_ID.value: run_id,
            db_handler.TableKeys.NAME.value: name,
            db_handler.TableKeys.VALUE.value: value,
        }
    )
    sqlite_insert.return_value.values.return_value.on_conflict_do_update.assert_not_called()
    execute.assert_called_once_with(statement)


def test_insert_with_conflic_mode_update(mocker, db_handler):
    table = mocker.ANY
    name = "test_name"
    value = "test_value"
    run_id = 1
    db_handler.run_id = run_id
    sqlite_insert = mocker.patch("qupled.database.sqlite_insert")
    statement = (
        sqlite_insert.return_value.values.return_value.on_conflict_do_update.return_value
    )
    execute = mocker.patch.object(db_handler, "_execute")
    db_handler._insert(table, name, value, DataBaseHandler.ConflictMode.UPDATE)
    sqlite_insert.assert_called_once_with(table)
    sqlite_insert.return_value.values.assert_called_once_with(
        {
            db_handler.TableKeys.RUN_ID.value: run_id,
            db_handler.TableKeys.NAME.value: name,
            db_handler.TableKeys.VALUE.value: value,
        }
    )
    sqlite_insert.return_value.values.return_value.on_conflict_do_update.assert_called_once_with(
        index_elements=[
            db_handler.TableKeys.RUN_ID.value,
            db_handler.TableKeys.NAME.value,
        ],
        set_={db_handler.TableKeys.VALUE.value: value},
    )
    execute.assert_called_once_with(statement)


def test_get(mocker, db_handler):
    run_id = 1
    names = ["a", "b"]
    table = mocker.MagicMock()
    table.c = mocker.MagicMock()
    select = mocker.patch("sqlalchemy.select")
    statement = select.return_value.where.return_value
    db_rows = [
        {db_handler.TableKeys.NAME.value: "a", db_handler.TableKeys.VALUE.value: 10},
        {db_handler.TableKeys.NAME.value: "b", db_handler.TableKeys.VALUE.value: 20},
    ]
    result = mocker.Mock()
    result.mappings.return_value.all.return_value = db_rows
    mocker.patch.object(db_handler, "_execute", return_value=result)
    sql_mapping = lambda x: x * 2
    actual = db_handler._get(table, run_id, names, sql_mapping)
    expected = {
        row[db_handler.TableKeys.NAME.value]: sql_mapping(
            row[db_handler.TableKeys.VALUE.value]
        )
        for row in db_rows
    }
    assert actual == expected
    select.assert_called_once_with(table)
    select.return_value.where.assert_called_once()
    db_handler._execute.assert_called_once_with(statement)


def test_execute(mocker, db_handler):
    statement = mocker.ANY
    connection = mocker.Mock()
    result = connection.execute.return_value
    engine = mocker.patch.object(db_handler.engine, "begin")
    engine.return_value.__enter__.return_value = connection
    result = db_handler._execute(statement)
    engine.assert_called_once()
    connection.execute.assert_called_once_with(statement)
    assert result == result


def test_to_bytes_with_float(db_handler):
    value = 3.14
    expected = struct.pack("d", value)
    assert db_handler._to_bytes(value) == expected


def test_to_bytes_with_numpy_array(db_handler):
    array = np.array([1.0, 2.0, 3.0])
    result = db_handler._to_bytes(array)
    assert isinstance(result, bytes)
    decompressed_result = blosc2.decompress(result)
    loaded = np.load(io.BytesIO(decompressed_result), allow_pickle=False)
    assert np.array_equal(loaded, array)


def test_to_bytes_with_invalid_type(db_handler):
    assert db_handler._to_bytes("invalid") is None


def test_from_bytes_with_float_bytes(db_handler):
    value = 42.0
    packed = struct.pack("d", value)
    assert db_handler._from_bytes(packed) == value


def test_from_bytes_with_numpy_bytes(db_handler):
    array = np.array([[1, 2], [3, 4]])
    arr_bytes = io.BytesIO()
    np.save(arr_bytes, array)
    compressed_arr_bytes = blosc2.compress(arr_bytes.getvalue())
    result = db_handler._from_bytes(compressed_arr_bytes)
    np.testing.assert_array_equal(result, array)


def test_from_bytes_with_invalid_bytes(db_handler):
    assert db_handler._from_bytes(b"not valid") is None


def test_to_json_valid_data(db_handler):
    data = {"a": 1, "b": 2}
    assert db_handler._to_json(data) == json.dumps(data)


def test_to_json_invalid_data(db_handler):
    class NotSerializable:
        pass

    assert db_handler._to_json(NotSerializable()) is None


def test_from_json_valid_data(db_handler):
    json_str = '{"x": 10, "y": 20}'
    assert db_handler._from_json(json_str) == json.loads(json_str)


def test_from_json_invalid_data(db_handler):
    assert db_handler._from_json("not a valid json") is None


# Functional tests


@pytest.fixture
def db_inputs():
    class Inputs:
        def __init__(self):
            self.theory = "theory"
            self.coupling = 0.5
            self.degeneracy = 2.0

    inputs = Inputs()
    yield inputs


@pytest.fixture
def db_results():
    class Results:
        def __init__(self):
            self.data = np.ndarray([1, 2, 3])

    results = Results()
    yield results


def test_insert_run_and_get_run_without_results(db_handler, db_inputs, db_results):
    db_handler.insert_run(db_inputs)
    run_data = db_handler.get_run(db_handler.run_id, None, None)
    run = run_data[DataBaseHandler.RUN_TABLE_NAME]
    inputs = run_data[DataBaseHandler.INPUT_TABLE_NAME]
    results = run_data[DataBaseHandler.RESULT_TABLE_NAME]
    assert run["theory"] == db_inputs.theory
    assert run["coupling"] == db_inputs.coupling
    assert run["degeneracy"] == db_inputs.degeneracy
    assert isinstance(run["date"], str)
    assert isinstance(run["time"], str)
    assert inputs["theory"] == db_inputs.theory
    assert inputs["coupling"] == db_inputs.coupling
    assert inputs["degeneracy"] == db_inputs.degeneracy
    assert results == {}


def test_insert_run_and_get_run_with_results(db_handler, db_inputs, db_results):
    db_handler.insert_run(db_inputs)
    db_handler.insert_results(db_results.__dict__)
    run_data = db_handler.get_run(db_handler.run_id, None, None)
    run = run_data[DataBaseHandler.RUN_TABLE_NAME]
    inputs = run_data[DataBaseHandler.INPUT_TABLE_NAME]
    results = run_data[DataBaseHandler.RESULT_TABLE_NAME]
    assert run["theory"] == db_inputs.theory
    assert run["coupling"] == db_inputs.coupling
    assert run["degeneracy"] == db_inputs.degeneracy
    assert isinstance(run["date"], str)
    assert isinstance(run["time"], str)
    assert run["status"] == DataBaseHandler.RunStatus.RUNNING.value
    assert inputs["theory"] == db_inputs.theory
    assert inputs["coupling"] == db_inputs.coupling
    assert inputs["degeneracy"] == db_inputs.degeneracy
    assert (results["data"] == db_results.data).all()


def test_insert_run_and_get_inputs(db_handler, db_inputs):
    db_handler.insert_run(db_inputs)
    inputs = db_handler.get_inputs(db_handler.run_id, None)
    assert inputs["theory"] == db_inputs.theory
    assert inputs["coupling"] == db_inputs.coupling
    assert inputs["degeneracy"] == db_inputs.degeneracy


def test_insert_run_and_get_results(db_handler, db_inputs, db_results):
    db_handler.insert_run(db_inputs)
    db_handler.insert_results(db_results.__dict__)
    results = db_handler.get_results(db_handler.run_id, None)
    assert (results["data"] == db_results.data).all()


def test_insert_inputs_without_run(db_handler, db_inputs):
    db_handler.insert_inputs(db_inputs.__dict__)
    inputs = db_handler.get_inputs(db_handler.run_id, None)
    assert inputs == {}


def test_insert_results_without_run(db_handler, db_results):
    db_handler.insert_results(db_results.__dict__)
    results = db_handler.get_results(db_handler.run_id, None)
    assert results == {}


def test_update_inputs_integrity_error(db_handler, db_inputs):
    db_handler.insert_run(db_inputs)
    with pytest.raises(sql.exc.IntegrityError):
        db_handler.insert_inputs(db_inputs.__dict__)


def test_update_results_default(db_handler, db_inputs, db_results):
    db_handler.insert_run(db_inputs)
    db_handler.insert_results(db_results.__dict__)
    with pytest.raises(sql.exc.IntegrityError):
        db_handler.insert_results(db_results.__dict__)


def test_update_results_allow_update(db_handler, db_inputs, db_results):
    db_handler.insert_run(db_inputs)
    db_handler.insert_results(db_results.__dict__)
    new_data = db_results.data + np.ones(db_results.data.shape)
    db_results.data = new_data
    db_handler.insert_results(
        db_results.__dict__, conflict_mode=DataBaseHandler.ConflictMode.UPDATE
    )
    results = db_handler.get_results(db_handler.run_id, None)
    assert (results["data"] == new_data).all()


def test_get_non_existing_run(db_handler, db_inputs):
    db_handler.insert_run(db_inputs)
    run_data = db_handler.get_run(db_handler.run_id + 1, None, None)
    assert run_data == {}


def test_get_non_existing_inputs(db_handler, db_inputs):
    db_handler.insert_run(db_inputs)
    inputs = db_handler.get_inputs(db_handler.run_id + 1, None)
    assert inputs == {}


def test_get_non_existing_results(db_handler, db_inputs):
    db_handler.insert_run(db_inputs)
    results = db_handler.get_results(db_handler.run_id + 1, None)
    assert results == {}


def test_insert_run_with_results_and_delete_run(db_handler, db_inputs, db_results):
    db_handler.insert_run(db_inputs)
    db_handler.insert_results(db_results.__dict__)
    run_data = db_handler.get_run(db_handler.run_id, None, None)
    assert run_data != {}
    db_handler.delete_run(db_handler.run_id)
    run_data = db_handler.get_run(db_handler.run_id, None, None)
    assert run_data == {}
