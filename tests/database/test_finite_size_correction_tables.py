import numpy as np
import pytest
import sqlalchemy as sql

from qupled.database.finite_size_correction_tables import (
    FiniteSizeCorrectionTables,
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
    tables = FiniteSizeCorrectionTables(engine)
    yield tables


@pytest.fixture
def tables_with_mock_build(mocker):
    mocker.patch("qupled.database.base_tables.BaseTables._build_tables")
    engine = sql.create_engine("sqlite:///:memory:")
    tables = FiniteSizeCorrectionTables(engine)
    yield tables


def test_inheritance():
    assert issubclass(FiniteSizeCorrectionTables, BaseTables)


def test_base_tables_initialization(tables):
    assert tables.engine is not None
    assert tables.run_table_name == RUN_TABLE_NAME
    assert tables.input_table_name == INPUT_TABLE_NAME
    assert tables.result_table_name == RESULT_TABLE_NAME
    assert tables.run_table is not None
    assert tables.input_table is not None
    assert tables.result_table is not None


def test_insert_inputs(mocker, tables):
    super_insert_inputs = mocker.patch(
        "qupled.database.base_tables.BaseTables.insert_inputs"
    )
    inputs = {"scheme": "scheme_name"}
    tables.insert_inputs(inputs)
    super_insert_inputs.assert_called_once_with({})


def test_build_run_table(mocker, tables_with_mock_build):
    create_table = mocker.patch.object(tables_with_mock_build, "_create_table")
    table = tables_with_mock_build._build_run_table()
    columns = {col.name for col in table.columns}
    expected_columns = {
        BaseTableKeys.PRIMARY_KEY.value,
        TableKeys.THEORY.value,
        TableKeys.COUPLING.value,
        TableKeys.DEGENERACY.value,
        TableKeys.NUMBER_OF_PARTICLES.value,
        TableKeys.SCHEME_RUN_ID.value,
        TableKeys.DATE.value,
        TableKeys.TIME.value,
        BaseTableKeys.STATUS.value,
    }
    create_table.assert_called_once_with(table)
    assert columns == expected_columns
    assert isinstance(table.c[BaseTableKeys.PRIMARY_KEY.value].type, sql.Integer)
    assert isinstance(table.c[TableKeys.THEORY.value].type, sql.String)
    assert isinstance(table.c[TableKeys.COUPLING.value].type, sql.Float)
    assert isinstance(table.c[TableKeys.DEGENERACY.value].type, sql.Float)
    assert isinstance(table.c[TableKeys.NUMBER_OF_PARTICLES.value].type, sql.Integer)
    # assert isinstance(table.c[TableKeys.SCHEME_RUN_ID.value].type, sql.Integer)
    assert isinstance(table.c[TableKeys.DATE.value].type, sql.String)
    assert isinstance(table.c[TableKeys.TIME.value].type, sql.String)
    assert isinstance(table.c[BaseTableKeys.STATUS.value].type, sql.String)
    assert not table.c[BaseTableKeys.PRIMARY_KEY.value].nullable
    assert not table.c[TableKeys.THEORY.value].nullable
    assert not table.c[TableKeys.COUPLING.value].nullable
    assert not table.c[TableKeys.DEGENERACY.value].nullable
    assert not table.c[TableKeys.NUMBER_OF_PARTICLES.value].nullable
    assert table.c[TableKeys.SCHEME_RUN_ID.value].nullable
    assert not table.c[TableKeys.DATE.value].nullable
    assert not table.c[TableKeys.TIME.value].nullable


def test_update_scheme_run_id_with_run_id(mocker, tables):
    tables.run_id = mocker.ANY
    scheme_run_id = mocker.ANY
    sql_update = mocker.patch("sqlalchemy.update")
    tables.run_table = mocker.MagicMock()
    statement = sql_update.return_value.where.return_value.values.return_value
    execute = mocker.patch.object(tables, "_execute")
    tables.update_scheme_run_id(scheme_run_id)
    sql_update.assert_called_once_with(tables.run_table)
    sql_update.return_value.where.return_value.values.assert_called_once_with(
        {TableKeys.SCHEME_RUN_ID.value: scheme_run_id}
    )
    execute.assert_called_once_with(statement)


def test_update_scheme_run_id_without_run_id(mocker, tables):
    tables.run_id = None
    status = 0
    sql_update = mocker.patch("sqlalchemy.update")
    execute = mocker.patch.object(tables, "_execute")
    tables.update_scheme_run_id(status)
    sql_update.assert_not_called()
    execute.assert_not_called()


def test_insert_run(mocker, tables):
    run_id = 1
    mock_datetime = mocker.patch(
        "qupled.database.finite_size_correction_tables.datetime"
    )
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
def fsc_inputs():
    class SchemeInputs:
        def __init__(self):
            self.theory = "theory"
            self.coupling = 1.0
            self.degeneracy = 1.0

    class Inputs:
        def __init__(self):
            self.number_of_particles = 1
            self.scheme = SchemeInputs()

    inputs = Inputs()
    yield inputs


@pytest.fixture
def fsc_results():
    class Results:
        def __init__(self):
            self.data = np.ndarray([1, 2, 3])

    results = Results()
    yield results


def test_insert_run_and_get_run_without_results(tables, fsc_inputs):
    tables.insert_run(fsc_inputs)
    run_data = tables.get_run(tables.run_id, None, None)
    run = run_data["run"]
    inputs = run_data["inputs"]
    results = run_data["results"]
    assert run["theory"] == fsc_inputs.scheme.theory
    assert run["coupling"] == fsc_inputs.scheme.coupling
    assert run["degeneracy"] == fsc_inputs.scheme.degeneracy
    assert run["number_of_particles"] == fsc_inputs.number_of_particles
    assert run["scheme_run_id"] is None
    assert isinstance(run["date"], str)
    assert isinstance(run["time"], str)
    assert inputs["number_of_particles"] == fsc_inputs.number_of_particles
    assert "scheme" not in inputs
    assert results == {}


def test_insert_run_and_get_run_with_results(tables, fsc_inputs, fsc_results):
    tables.insert_run(fsc_inputs)
    tables.insert_results(fsc_results.__dict__)
    run_data = tables.get_run(tables.run_id, None, None)
    run = run_data["run"]
    inputs = run_data["inputs"]
    results = run_data["results"]
    assert run["theory"] == fsc_inputs.scheme.theory
    assert run["coupling"] == fsc_inputs.scheme.coupling
    assert run["degeneracy"] == fsc_inputs.scheme.degeneracy
    assert run["number_of_particles"] == fsc_inputs.number_of_particles
    assert run["scheme_run_id"] is None
    assert isinstance(run["date"], str)
    assert isinstance(run["time"], str)
    assert run["status"] == RunStatus.RUNNING.value
    assert inputs["number_of_particles"] == fsc_inputs.number_of_particles
    assert "scheme" not in inputs
    assert np.array_equal(results["data"], fsc_results.data, equal_nan=True)


def test_insert_run_and_get_inputs(tables, fsc_inputs):
    tables.insert_run(fsc_inputs)
    inputs = tables.get_inputs(tables.run_id, None)
    assert inputs["number_of_particles"] == fsc_inputs.number_of_particles


def test_insert_run_and_get_results(tables, fsc_inputs, fsc_results):
    tables.insert_run(fsc_inputs)
    tables.insert_results(fsc_results.__dict__)
    results = tables.get_results(tables.run_id, None)
    assert np.array_equal(results["data"], fsc_results.data, equal_nan=True)


def test_insert_inputs_without_run(tables, fsc_inputs):
    tables.insert_inputs(fsc_inputs.__dict__)
    inputs = tables.get_inputs(tables.run_id, None)
    assert inputs == {}


def test_insert_results_without_run(tables, fsc_results):
    tables.insert_results(fsc_results.__dict__)
    results = tables.get_results(tables.run_id, None)
    assert results == {}


def test_update_inputs_integrity_error(tables, fsc_inputs):
    tables.insert_run(fsc_inputs)
    with pytest.raises(sql.exc.IntegrityError):
        tables.insert_inputs(fsc_inputs.__dict__)


def test_update_results_default(tables, fsc_inputs, fsc_results):
    tables.insert_run(fsc_inputs)
    tables.insert_results(fsc_results.__dict__)
    with pytest.raises(sql.exc.IntegrityError):
        tables.insert_results(fsc_results.__dict__)


def test_update_results_allow_update(tables, fsc_inputs, fsc_results):
    tables.insert_run(fsc_inputs)
    tables.insert_results(fsc_results.__dict__)
    new_data = fsc_results.data + np.ones(fsc_results.data.shape)
    fsc_results.data = new_data
    tables.insert_results(fsc_results.__dict__, conflict_mode=ConflictMode.UPDATE)
    results = tables.get_results(tables.run_id, None)
    assert np.array_equal(results["data"], new_data, equal_nan=True)


def test_get_non_existing_run(tables, fsc_inputs):
    tables.insert_run(fsc_inputs)
    run_data = tables.get_run(tables.run_id + 1, None, None)
    assert run_data == {}


def test_get_non_existing_inputs(tables, fsc_inputs):
    tables.insert_run(fsc_inputs)
    inputs = tables.get_inputs(tables.run_id + 1, None)
    assert inputs == {}


def test_get_non_existing_results(tables, fsc_inputs):
    tables.insert_run(fsc_inputs)
    results = tables.get_results(tables.run_id + 1, None)
    assert results == {}


def test_insert_run_with_results_and_delete_run(tables, fsc_inputs, fsc_results):
    tables.insert_run(fsc_inputs)
    tables.insert_results(fsc_results.__dict__)
    run_data = tables.get_run(tables.run_id, None, None)
    assert run_data != {}
    tables.delete_run(tables.run_id)
    run_data = tables.get_run(tables.run_id, None, None)
    assert run_data == {}
