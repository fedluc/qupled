import numpy as np
import pytest
import sqlalchemy as sql

from qupled.database.scheme_tables import (
    SchemeTables,
    BaseTables,
    BaseTableKeys,
    RunStatus,
    TableKeys,
    INPUT_TABLE_NAME,
    RESULT_TABLE_NAME,
    RUN_TABLE_NAME,
)

from qupled.database.base_tables import ConflictMode

# Unit tests

@pytest.fixture
def tables():
    engine = sql.create_engine("sqlite:///:memory:")
    tables = SchemeTables(engine)
    yield tables

@pytest.fixture
def tables_with_mock_build(mocker):
    mocker.patch("qupled.database.base_tables.BaseTables._build_tables")
    engine = sql.create_engine("sqlite:///:memory:")
    tables = SchemeTables(engine)
    yield tables


def test_inheritance():
    assert issubclass(SchemeTables, BaseTables)


def test_base_tables_initialization(tables):
    assert tables.engine is not None
    assert tables.run_table_name == RUN_TABLE_NAME
    assert tables.input_table_name == INPUT_TABLE_NAME
    assert tables.result_table_name == RESULT_TABLE_NAME
    assert tables.run_table is not None
    assert tables.input_table is not None
    assert tables.result_table is not None

def test_build_run_table(mocker, tables_with_mock_build):
    create_table = mocker.patch.object(tables_with_mock_build, "_create_table")
    table = tables_with_mock_build._build_run_table()
    columns = {col.name for col in table.columns}
    expected_columns = {
        BaseTableKeys.PRIMARY_KEY.value,
        TableKeys.THEORY.value,
        TableKeys.COUPLING.value,
        TableKeys.DEGENERACY.value,
        TableKeys.DATE.value,
        TableKeys.TIME.value,
        BaseTableKeys.STATUS.value
    }
    create_table.assert_called_once_with(table)
    assert columns == expected_columns
    assert isinstance(table.c[BaseTableKeys.PRIMARY_KEY.value].type, sql.Integer)
    assert isinstance(table.c[TableKeys.THEORY.value].type, sql.String)
    assert isinstance(table.c[TableKeys.COUPLING.value].type, sql.Float)
    assert isinstance(table.c[TableKeys.DEGENERACY.value].type, sql.Float)
    assert isinstance(table.c[TableKeys.DATE.value].type, sql.String)
    assert isinstance(table.c[TableKeys.TIME.value].type, sql.String)
    assert isinstance(table.c[BaseTableKeys.STATUS.value].type, sql.String)
    assert not table.c[BaseTableKeys.PRIMARY_KEY.value].nullable
    assert not table.c[TableKeys.THEORY.value].nullable
    assert not table.c[TableKeys.COUPLING.value].nullable
    assert not table.c[TableKeys.DEGENERACY.value].nullable
    assert not table.c[TableKeys.DATE.value].nullable
    assert not table.c[TableKeys.TIME.value].nullable

def test_delete_run(mocker, tables):
    run_id = mocker.ANY
    delete_blob_data_on_disk = mocker.patch("qupled.database.scheme_tables.delete_blob_data_on_disk")
    super_delete_run =  mocker.patch("qupled.database.base_tables.BaseTables.delete_run")
    tables.delete_run(run_id)
    delete_blob_data_on_disk.assert_called_once_with(tables.engine.url.database, run_id)
    super_delete_run.assert_called_once_with(run_id)
    
def test_insert_run(mocker, tables):
    run_id = 1
    mock_datetime = mocker.patch("qupled.database.scheme_tables.datetime")
    mock_datetime.now.return_value = mocker.Mock()
    inputs = mocker.Mock()
    status = mocker.Mock()
    data = mocker.ANY
    sql_insert = mocker.patch("sqlalchemy.insert")
    statement = sql_insert.return_value.values.return_value
    result = mocker.Mock()
    execute = mocker.patch.object(tables, "_execute", return_value=result)
    result.inserted_primary_key = [run_id]
    tables._insert_run(inputs, status)
    mock_datetime.now.assert_called_once()
    sql_insert.assert_called_once_with(tables.run_table)
    sql_insert.return_value.values.assert_called_once_with(data)
    execute.assert_called_once_with(statement)
    assert tables.run_id == run_id


# Functional tests

@pytest.fixture
def scheme_inputs():
    class Inputs:
        def __init__(self):
            self.theory = "theory"
            self.coupling = 1.0
            self.degeneracy = 1.0

    inputs = Inputs()
    yield inputs


@pytest.fixture
def scheme_results():
    class Results:
        def __init__(self):
            self.data = np.ndarray([1, 2, 3])

    results = Results()
    yield results


def test_insert_run_and_get_run_without_results(tables, scheme_inputs):
    tables.insert_run(scheme_inputs)
    run_data = tables.get_run(tables.run_id, None, None)
    run = run_data[RUN_TABLE_NAME]
    inputs = run_data[INPUT_TABLE_NAME]
    results = run_data[RESULT_TABLE_NAME]
    assert run["theory"] == scheme_inputs.theory
    assert run["coupling"] == scheme_inputs.coupling
    assert run["degeneracy"] == scheme_inputs.degeneracy
    assert isinstance(run["date"], str)
    assert isinstance(run["time"], str)
    assert inputs["theory"] == scheme_inputs.theory
    assert inputs["coupling"] == scheme_inputs.coupling
    assert inputs["degeneracy"] == scheme_inputs.degeneracy
    assert results == {}


def test_insert_run_and_get_run_with_results(tables, scheme_inputs, scheme_results):
    tables.insert_run(scheme_inputs)
    tables.insert_results(scheme_results.__dict__)
    run_data = tables.get_run(tables.run_id, None, None)
    run = run_data[RUN_TABLE_NAME]
    inputs = run_data[INPUT_TABLE_NAME]
    results = run_data[RESULT_TABLE_NAME]
    assert run["theory"] == scheme_inputs.theory
    assert run["coupling"] == scheme_inputs.coupling
    assert run["degeneracy"] == scheme_inputs.degeneracy
    assert isinstance(run["date"], str)
    assert isinstance(run["time"], str)
    assert run["status"] == RunStatus.RUNNING.value
    assert inputs["theory"] == scheme_inputs.theory
    assert inputs["coupling"] == scheme_inputs.coupling
    assert inputs["degeneracy"] == scheme_inputs.degeneracy
    assert (results["data"] == scheme_results.data).all()


def test_insert_run_and_get_inputs(tables, scheme_inputs):
    tables.insert_run(scheme_inputs)
    inputs = tables.get_inputs(tables.run_id, None)
    assert inputs["theory"] == scheme_inputs.theory
    assert inputs["coupling"] == scheme_inputs.coupling
    assert inputs["degeneracy"] == scheme_inputs.degeneracy


def test_insert_run_and_get_results(tables, scheme_inputs, scheme_results):
    tables.insert_run(scheme_inputs)
    tables.insert_results(scheme_results.__dict__)
    results = tables.get_results(tables.run_id, None)
    assert (results["data"] == scheme_results.data).all()


def test_insert_inputs_without_run(tables, scheme_inputs):
    tables.insert_inputs(scheme_inputs.__dict__)
    inputs = tables.get_inputs(tables.run_id, None)
    assert inputs == {}


def test_insert_results_without_run(tables, scheme_results):
    tables.insert_results(scheme_results.__dict__)
    results = tables.get_results(tables.run_id, None)
    assert results == {}


def test_update_inputs_integrity_error(tables, scheme_inputs):
    tables.insert_run(scheme_inputs)
    with pytest.raises(sql.exc.IntegrityError):
        tables.insert_inputs(scheme_inputs.__dict__)


def test_update_results_default(tables, scheme_inputs, scheme_results):
    tables.insert_run(scheme_inputs)
    tables.insert_results(scheme_results.__dict__)
    with pytest.raises(sql.exc.IntegrityError):
        tables.insert_results(scheme_results.__dict__)


def test_update_results_allow_update(tables, scheme_inputs, scheme_results):
    tables.insert_run(scheme_inputs)
    tables.insert_results(scheme_results.__dict__)
    new_data = scheme_results.data + np.ones(scheme_results.data.shape)
    scheme_results.data = new_data
    tables.insert_results(
        scheme_results.__dict__, conflict_mode=ConflictMode.UPDATE
    )
    results = tables.get_results(tables.run_id, None)
    assert (results["data"] == new_data).all()


def test_get_non_existing_run(tables, scheme_inputs):
    tables.insert_run(scheme_inputs)
    run_data = tables.get_run(tables.run_id + 1, None, None)
    assert run_data == {}


def test_get_non_existing_inputs(tables, scheme_inputs):
    tables.insert_run(scheme_inputs)
    inputs = tables.get_inputs(tables.run_id + 1, None)
    assert inputs == {}


def test_get_non_existing_results(tables, scheme_inputs):
    tables.insert_run(scheme_inputs)
    results = tables.get_results(tables.run_id + 1, None)
    assert results == {}


def test_insert_run_with_results_and_delete_run(tables, scheme_inputs, scheme_results):
    tables.insert_run(scheme_inputs)
    tables.insert_results(scheme_results.__dict__)
    run_data = tables.get_run(tables.run_id, None, None)
    assert run_data != {}
    tables.delete_run(tables.run_id)
    run_data = tables.get_run(tables.run_id, None, None)
    assert run_data == {}
