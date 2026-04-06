import io
import json
import struct

import blosc2
import numpy as np
import pytest
import sqlalchemy as sql

from qupled.database.base_tables import BaseTables, ConflictMode, RunStatus, TableKeys


@pytest.fixture
def tables():
    engine = sql.create_engine("sqlite:///:memory:")
    tables = BaseTables(engine, "run_table", "input_table", "result_table")
    yield tables


@pytest.mark.unit
def test_base_tables_initialization(tables):
    assert tables.engine is not None
    assert tables.run_table_name == "run_table"
    assert tables.input_table_name == "input_table"
    assert tables.result_table_name == "result_table"
    assert tables.table_metadata is not None
    assert tables.run_table is None
    assert tables.input_table is None
    assert tables.result_table is None
    assert tables.run_id is None


@pytest.mark.unit
def test_insert_run(mocker, tables):
    insert_run = mocker.patch.object(tables, "_insert_run")
    insert_inputs = mocker.patch.object(tables, "insert_inputs")
    inputs = mocker.ANY
    tables.insert_run(inputs)
    insert_run.assert_called_once_with(inputs, RunStatus.RUNNING)
    insert_inputs.assert_called_once_with(inputs.__dict__)


@pytest.mark.unit
def test_insert_inputs(mocker, tables):
    tables.run_id = mocker.ANY
    insert_data_from_dict = mocker.patch.object(tables, "_insert_data_from_dict")
    mocker.patch.object(tables, "_to_json", side_effect=lambda x: f"json({x})")
    inputs = mocker.ANY
    tables.insert_inputs(inputs)
    insert_data_from_dict.assert_called_once()
    called_table, called_inputs, called_mapper = insert_data_from_dict.call_args[0]
    assert called_table == tables.input_table
    assert called_inputs == inputs
    assert callable(called_mapper)
    assert called_mapper("value1") == "json(value1)"
    assert called_mapper(123) == "json(123)"


@pytest.mark.unit
def test_insert_inputs_without_run_id(mocker, tables):
    insert_from_dict = mocker.patch.object(tables, "_insert_data_from_dict")
    mocker.patch.object(tables, "_to_json", side_effect=lambda x: f"json({x})")
    tables.insert_inputs(mocker.ANY)
    insert_from_dict.assert_not_called()


@pytest.mark.unit
def test_insert_results(mocker, tables):
    tables.run_id = mocker.ANY
    insert_from_dict = mocker.patch.object(tables, "_insert_data_from_dict")
    mocker.patch.object(tables, "_to_bytes", side_effect=lambda x: f"bytes({x})")
    results = mocker.ANY
    tables.insert_results(results)
    insert_from_dict.assert_called_once()
    called_table, called_results, called_mapper, called_conflict_mode = (
        insert_from_dict.call_args[0]
    )
    assert called_table == tables.result_table
    assert called_results == results
    assert callable(called_mapper)
    assert called_mapper("value1") == "bytes(value1)"
    assert called_mapper(123) == "bytes(123)"
    assert called_conflict_mode == ConflictMode.FAIL


@pytest.mark.unit
def test_insert_results_without_run_id(mocker, tables):
    insert_from_dict = mocker.patch.object(tables, "_insert_data_from_dict")
    mocker.patch.object(tables, "_to_json", side_effect=lambda x: f"json({x})")
    tables.insert_results(mocker.ANY)
    insert_from_dict.assert_not_called()


@pytest.mark.unit
def test_inspect_runs(mocker, tables):
    sql_select = mocker.patch("sqlalchemy.select")
    mock_runs = [{"id": 1, "status": "done"}, {"id": 2, "status": "failed"}]
    execute = mocker.patch.object(tables, "_execute")
    mock_result = mocker.Mock()
    mock_result.mappings.return_value.all.return_value = mock_runs
    execute.return_value = mock_result
    runs = tables.inspect_runs()
    sql_select.assert_called_once_with(tables.run_table)
    mock_result.mappings.return_value.all.assert_called_once()
    assert runs == mock_runs


@pytest.mark.unit
def test_update_run_status_with_run_id(mocker, tables):
    tables.run_id = 1
    status = RunStatus.SUCCESS
    sql_update = mocker.patch("sqlalchemy.update")
    tables.run_table = mocker.MagicMock()
    statement = sql_update.return_value.where.return_value.values.return_value
    execute = mocker.patch.object(tables, "_execute")
    tables.update_run_status(status)
    sql_update.assert_called_once_with(tables.run_table)
    sql_update.return_value.where.return_value.values.assert_called_once_with(
        {TableKeys.STATUS.value: status.value}
    )
    execute.assert_called_once_with(statement)


@pytest.mark.unit
def test_update_run_status_without_run_id(mocker, tables):
    tables.run_id = None
    status = 0
    sql_update = mocker.patch("sqlalchemy.update")
    execute = mocker.patch.object(tables, "_execute")
    tables.update_run_status(status)
    sql_update.assert_not_called()
    execute.assert_not_called()


@pytest.mark.unit
def test_get_run_with_existing_run(mocker, tables):
    run_id = 1
    sql_select = mocker.patch("sqlalchemy.select")
    execute = mocker.patch.object(tables, "_execute")
    statement = sql_select.return_value.where.return_value
    tables.run_table = mocker.MagicMock()
    mock_result = mocker.Mock()
    mock_result.mappings.return_value.first.return_value = {"key": "value"}
    execute.return_value = mock_result
    inputs = mocker.ANY
    get_inputs = mocker.patch.object(tables, "get_inputs", return_value=inputs)
    results = mocker.ANY
    get_results = mocker.patch.object(tables, "get_results", return_value=results)
    run = tables.get_run(run_id, None, None)
    sql_select.assert_called_once_with(tables.run_table)
    get_inputs.assert_called_once_with(run_id, names=None)
    get_results.assert_called_once_with(run_id, names=None)
    execute.assert_called_once_with(statement)
    mock_result.mappings.return_value.first.assert_called_once()
    assert run.run == {"key": "value"}
    assert run.inputs == inputs
    assert run.results == results


@pytest.mark.unit
def test_get_run_with_non_existing_run(mocker, tables):
    run_id = 1
    sql_select = mocker.patch("sqlalchemy.select")
    statement = sql_select.return_value.where.return_value
    tables.run_table = mocker.MagicMock()
    execute = mocker.patch.object(tables, "_execute")
    mock_result = mocker.Mock()
    mock_result.mappings.return_value.first.return_value = None
    execute.return_value = mock_result
    run = tables.get_run(run_id, None, None)
    sql_select.assert_called_once_with(tables.run_table)
    execute.assert_called_once_with(statement)
    mock_result.mappings.return_value.first.assert_called_once()
    assert run is None


@pytest.mark.unit
def test_get_inputs(mocker, tables):
    run_id = 1
    names = ["name"]
    expected_inputs = {"key": "value"}
    get = mocker.patch.object(tables, "_get_data", return_value=expected_inputs)
    mocker.patch.object(tables, "_from_json", side_effect=lambda x: f"from_json({x})")
    inputs = tables.get_inputs(run_id, names)
    get.assert_called_once()
    called_table, called_run_id, called_names, called_mapper = get.call_args[0]
    assert called_table == tables.input_table
    assert called_run_id == run_id
    assert called_names == names
    assert callable(called_mapper)
    assert called_mapper("value1") == "from_json(value1)"
    assert called_mapper(123) == "from_json(123)"
    assert inputs == expected_inputs


@pytest.mark.unit
def test_get_results(mocker, tables):
    run_id = 1
    names = ["name"]
    expected_results = {"key": "value"}
    get = mocker.patch.object(tables, "_get_data", return_value=expected_results)
    mocker.patch.object(tables, "_from_bytes", side_effect=lambda x: f"from_bytes({x})")
    results = tables.get_results(run_id, names)
    get.assert_called_once()
    called_table, called_run_id, called_names, called_mapper = get.call_args[0]
    assert called_table == tables.result_table
    assert called_run_id == run_id
    assert called_names == names
    assert callable(called_mapper)
    assert called_mapper("value1") == "from_bytes(value1)"
    assert called_mapper(123) == "from_bytes(123)"
    assert results == expected_results


@pytest.mark.unit
def test_delete_run(mocker, tables):
    run_id = 1
    tables.run_table = mocker.MagicMock()
    tables.run_table.result_value.c = mocker.ANY
    sql_delete = mocker.patch("sqlalchemy.delete")
    statement = sql_delete.return_value.where.return_value
    execute = mocker.patch.object(tables, "_execute")
    tables.delete_run(run_id)
    sql_delete.assert_called_once_with(tables.run_table)
    execute.assert_called_once_with(statement)


@pytest.mark.unit
def test_build_tables(mocker, tables):
    mock_run_table = mocker.Mock()
    mock_input_table = mocker.Mock()
    mock_result_table = mocker.Mock()
    build_run_table = mocker.patch.object(
        tables, "_build_run_table", return_value=mock_run_table
    )
    build_inputs_table = mocker.patch.object(
        tables, "_build_inputs_table", return_value=mock_input_table
    )
    build_results_table = mocker.patch.object(
        tables, "_build_results_table", return_value=mock_result_table
    )
    tables._build_tables()
    build_run_table.assert_called_once()
    build_inputs_table.assert_called_once()
    build_results_table.assert_called_once()
    assert build_run_table.return_value == mock_run_table
    assert build_inputs_table.return_value == mock_input_table
    assert build_results_table.return_value == mock_result_table
    assert tables.run_table == mock_run_table
    assert tables.input_table == mock_input_table
    assert tables.result_table == mock_result_table


@pytest.mark.unit
def test_build_run_table(tables):
    with pytest.raises(NotImplementedError) as excinfo:
        tables._build_run_table()
    assert excinfo.value.args[0] == "This method should be implemented in a subclass."


@pytest.mark.unit
def test_build_inputs_table(mocker, tables):
    mock_table = mocker.Mock()
    build_data_table = mocker.patch.object(
        tables, "_build_data_table", return_value=mock_table
    )
    tables._build_inputs_table()
    build_data_table.assert_called_once_with(tables.input_table_name, sql.JSON)
    assert build_data_table.return_value == mock_table


@pytest.mark.unit
def test_build_results_table(mocker, tables):
    mock_table = mocker.Mock()
    build_data_table = mocker.patch.object(
        tables, "_build_data_table", return_value=mock_table
    )
    tables._build_results_table()
    build_data_table.assert_called_once_with(tables.result_table_name, sql.LargeBinary)
    assert build_data_table.return_value == mock_table


@pytest.mark.unit
def test_build_data_table(mocker, tables):
    tables.run_table = sql.Table(
        tables.run_table_name,
        tables.table_metadata,
        sql.Column(
            TableKeys.PRIMARY_KEY.value,
            sql.Integer,
            primary_key=True,
            autoincrement=True,
        ),
    )
    create_table = mocker.patch.object(tables, "_create_table")
    table = tables._build_data_table("table_name", sql.JSON)
    columns = {col.name for col in table.columns}
    indexes = {index.name for index in table.indexes}
    expected_columns = {
        TableKeys.RUN_ID.value,
        TableKeys.NAME.value,
        TableKeys.VALUE.value,
    }
    expected_indexes = {"idx_table_name_run_id", "idx_table_name_name"}
    create_table.assert_called_once_with(table)
    assert columns == expected_columns
    assert indexes == expected_indexes
    assert isinstance(table.c[TableKeys.RUN_ID.value].type, sql.Integer)
    assert isinstance(table.c[TableKeys.NAME.value].type, sql.String)
    assert isinstance(table.c[TableKeys.VALUE.value].type, sql.JSON)
    assert not table.c[TableKeys.RUN_ID.value].nullable
    assert not table.c[TableKeys.NAME.value].nullable
    assert table.c[TableKeys.VALUE.value].nullable


@pytest.mark.unit
def test_create_table(mocker, tables):
    mock_table = mocker.Mock()
    create = mocker.patch.object(mock_table, "create")
    tables._create_table(mock_table)
    create.assert_called_once_with(tables.engine, checkfirst=True)


@pytest.mark.unit
def test_insert_run(mocker, tables):
    with pytest.raises(NotImplementedError) as excinfo:
        tables._insert_run(mocker.ANY, mocker.ANY)
    assert excinfo.value.args[0] == "This method should be implemented in a subclass."


@pytest.mark.unit
def test_insert_data_from_dict_with_valid_data(mocker, tables):
    table = mocker.ANY
    data = {"key1": "value1", "key2": "value2"}
    sql_mapping = mocker.Mock(side_effect=lambda x: f"mapped({x})")
    insert_data = mocker.patch.object(tables, "_insert_data")
    tables._insert_data_from_dict(table, data, sql_mapping)
    sql_mapping.assert_has_calls([mocker.call("value1"), mocker.call("value2")])
    insert_data.assert_has_calls(
        [
            mocker.call(table, "key1", "mapped(value1)", mocker.ANY),
            mocker.call(table, "key2", "mapped(value2)", mocker.ANY),
        ]
    )
    assert insert_data.call_count == 2


@pytest.mark.unit
def test_insert_data_from_dict_with_empty_data(mocker, tables):
    data = {}
    sql_mapping = mocker.Mock()
    insert_data = mocker.patch.object(tables, "_insert_data")
    tables._insert_data_from_dict(mocker.ANY, data, sql_mapping)
    sql_mapping.assert_not_called()
    insert_data.assert_not_called()


@pytest.mark.unit
def test_insert_data_with_conflict_mode_fail(mocker, tables):
    table = mocker.ANY
    name = "test_name"
    value = "test_value"
    run_id = 1
    tables.run_id = run_id
    sqlite_insert = mocker.patch("qupled.database.base_tables.sqlite_insert")
    statement = sqlite_insert.return_value.values.return_value
    execute = mocker.patch.object(tables, "_execute")
    tables._insert_data(table, name, value)
    sqlite_insert.assert_called_once_with(table)
    sqlite_insert.return_value.values.assert_called_once_with(
        {
            TableKeys.RUN_ID.value: run_id,
            TableKeys.NAME.value: name,
            TableKeys.VALUE.value: value,
        }
    )
    sqlite_insert.return_value.values.return_value.on_conflict_do_update.assert_not_called()
    execute.assert_called_once_with(statement)


@pytest.mark.unit
def test_insert_data_with_conflict_mode_update(mocker, tables):
    table = mocker.ANY
    name = "test_name"
    value = "test_value"
    run_id = 1
    tables.run_id = run_id
    sqlite_insert = mocker.patch("qupled.database.base_tables.sqlite_insert")
    statement = (
        sqlite_insert.return_value.values.return_value.on_conflict_do_update.return_value
    )
    execute = mocker.patch.object(tables, "_execute")
    tables._insert_data(table, name, value, ConflictMode.UPDATE)
    sqlite_insert.assert_called_once_with(table)
    sqlite_insert.return_value.values.assert_called_once_with(
        {
            TableKeys.RUN_ID.value: run_id,
            TableKeys.NAME.value: name,
            TableKeys.VALUE.value: value,
        }
    )
    sqlite_insert.return_value.values.return_value.on_conflict_do_update.assert_called_once_with(
        index_elements=[
            TableKeys.RUN_ID.value,
            TableKeys.NAME.value,
        ],
        set_={TableKeys.VALUE.value: value},
    )
    execute.assert_called_once_with(statement)


@pytest.mark.unit
def test_get_data(mocker, tables):
    run_id = 1
    names = ["a", "b"]
    table = mocker.MagicMock()
    table.c = mocker.MagicMock()
    select = mocker.patch("sqlalchemy.select")
    statement = select.return_value.where.return_value
    db_rows = [
        {TableKeys.NAME.value: "a", TableKeys.VALUE.value: 10},
        {TableKeys.NAME.value: "b", TableKeys.VALUE.value: 20},
    ]
    result = mocker.Mock()
    result.mappings.return_value.all.return_value = db_rows
    mocker.patch.object(tables, "_execute", return_value=result)
    sql_mapping = lambda x: x * 2
    actual = tables._get_data(table, run_id, names, sql_mapping)
    expected = {
        row[TableKeys.NAME.value]: sql_mapping(row[TableKeys.VALUE.value])
        for row in db_rows
    }
    assert actual == expected
    select.assert_called_once_with(table)
    select.return_value.where.assert_called_once()
    tables._execute.assert_called_once_with(statement)


@pytest.mark.unit
def test_execute(mocker, tables):
    statement = mocker.ANY
    connection = mocker.Mock()
    result = connection.execute.return_value
    engine = mocker.patch.object(tables.engine, "begin")
    engine.return_value.__enter__.return_value = connection
    result = tables._execute(statement)
    engine.assert_called_once()
    connection.execute.assert_called_once_with(statement)
    assert result == result


@pytest.mark.unit
def test_to_bytes_with_float(tables):
    value = 3.14
    expected = struct.pack("d", value)
    assert tables._to_bytes(value) == expected


@pytest.mark.unit
def test_to_bytes_with_numpy_array(tables):
    array = np.array([1.0, 2.0, 3.0])
    result = tables._to_bytes(array)
    assert isinstance(result, bytes)
    decompressed_result = blosc2.decompress(result)
    loaded = np.load(io.BytesIO(decompressed_result), allow_pickle=False)
    assert np.array_equal(loaded, array)


@pytest.mark.unit
def test_to_bytes_with_invalid_type(tables):
    assert tables._to_bytes("invalid") is None


@pytest.mark.unit
def test_from_bytes_with_float_bytes(tables):
    value = 42.0
    packed = struct.pack("d", value)
    assert tables._from_bytes(packed) == value


@pytest.mark.unit
def test_from_bytes_with_numpy_bytes(tables):
    array = np.array([[1, 2], [3, 4]])
    arr_bytes = io.BytesIO()
    np.save(arr_bytes, array)
    compressed_arr_bytes = blosc2.compress(arr_bytes.getvalue())
    result = tables._from_bytes(compressed_arr_bytes)
    np.testing.assert_array_equal(result, array)


@pytest.mark.unit
def test_from_bytes_with_invalid_bytes(tables):
    assert tables._from_bytes(b"not valid") is None


@pytest.mark.unit
def test_to_json_valid_data(tables):
    data = {"a": 1, "b": 2}
    assert tables._to_json(data) == json.dumps(data)


@pytest.mark.unit
def test_to_json_invalid_data(tables):
    class NotSerializable:
        pass

    assert tables._to_json(NotSerializable()) is None


@pytest.mark.unit
def test_from_json_valid_data(tables):
    json_str = '{"x": 10, "y": 20}'
    assert tables._from_json(json_str) == json.loads(json_str)


@pytest.mark.unit
def test_from_json_invalid_data(tables):
    assert tables._from_json("not a valid json") is None
