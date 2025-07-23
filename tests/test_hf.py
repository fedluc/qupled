import pytest
import numpy as np
from unittest.mock import PropertyMock

from qupled.database import DataBaseHandler
import qupled.hf as hf
import qupled.native as native


@pytest.fixture
def inputs():
    return hf.Input(coupling=1.0, degeneracy=2.0)


@pytest.fixture
def results():
    return hf.Result()


@pytest.fixture
def scheme(mocker):
    scheme = hf.Solver()
    scheme.db_handler = mocker.Mock()
    return scheme


def test_native_to_run_status_mapping():
    assert hf.Solver.NATIVE_TO_RUN_STATUS[0] == DataBaseHandler.RunStatus.SUCCESS
    assert hf.Solver.NATIVE_TO_RUN_STATUS[1] == DataBaseHandler.RunStatus.FAILED
    assert len(hf.Solver.NATIVE_TO_RUN_STATUS) == 2


def test_hf_initialization():
    scheme = hf.Solver()
    assert scheme.inputs is None
    assert isinstance(scheme.results, hf.Result)
    assert isinstance(scheme.db_handler, DataBaseHandler)
    assert scheme.native_scheme_cls == native.HF
    assert scheme.native_inputs_cls, native.Input
    assert scheme.native_scheme_status is None


def test_run_id(scheme):
    run_id = "run_id"
    scheme.db_handler.run_id = run_id
    assert scheme.run_id == run_id


def test_compute(scheme, inputs, mocker):
    add_run_to_database = mocker.patch.object(scheme, "_add_run_to_database")
    compute_native = mocker.patch.object(scheme, "_compute_native")
    save = mocker.patch.object(scheme, "_save")
    scheme.compute(inputs)
    assert scheme.inputs is not None
    add_run_to_database.assert_called_once()
    compute_native.assert_called_once()
    save.assert_called_once()


def test_add_run_to_database(scheme, mocker):
    mocker.patch.object(hf.Solver, "run_id", new_callable=PropertyMock).return_value = (
        "mocked-run-id"
    )
    scheme.inputs = mocker.Mock()
    scheme._add_run_to_database()
    scheme.db_handler.insert_run.assert_called_once_with(scheme.inputs)
    assert scheme.inputs.database_info.run_id == scheme.run_id


def test_compute_native(scheme, inputs, mocker):
    native_input = mocker.Mock()
    native_inputs_cls = mocker.patch.object(
        scheme, "native_inputs_cls", return_value=native_input
    )
    to_native = mocker.patch("qupled.hf.Input.to_native")
    native_scheme = mocker.Mock()
    native_scheme_cls = mocker.patch.object(
        scheme, "native_scheme_cls", return_value=native_scheme
    )
    from_native = mocker.patch("qupled.hf.Result.from_native")
    native_scheme.compute.return_value = "mocked-status"
    scheme.inputs = inputs
    scheme._compute_native()
    native_inputs_cls.assert_called_once()
    to_native.assert_called_once_with(native_input)
    native_scheme_cls.assert_called_once_with(native_input)
    native_scheme.compute.assert_called_once()
    assert scheme.native_scheme_status == "mocked-status"
    from_native.assert_called_once_with(native_scheme)


def test_save(scheme, results, mocker):
    scheme.results = results
    scheme.native_scheme_status = mocker.Mock()
    scheme._save()
    scheme.db_handler.update_run_status.assert_called_once_with(
        hf.Solver.NATIVE_TO_RUN_STATUS.get(
            scheme.native_scheme_status, DataBaseHandler.RunStatus.FAILED
        )
    )
    scheme.db_handler.insert_results.assert_called_once_with(scheme.results.__dict__)


def test_compute_rdf_with_default_grid(scheme, inputs, results, mocker):
    compute_rdf = mocker.patch("qupled.hf.Result.compute_rdf")
    scheme.results = results
    scheme.inputs = inputs
    scheme.compute_rdf()
    compute_rdf.assert_called_once_with(scheme.inputs.dimension, None)
    scheme.db_handler.insert_results.assert_called_once_with(
        {
            "rdf": scheme.results.rdf,
            "rdf_grid": scheme.results.rdf_grid,
        },
        conflict_mode=DataBaseHandler.ConflictMode.UPDATE,
    )


def test_compute_rdf_with_custom_grid(scheme, inputs, results, mocker):
    compute_rdf = mocker.patch("qupled.hf.Result.compute_rdf")
    scheme.results = results
    scheme.inputs = inputs
    rdf_grid = np.array([1, 2, 3])
    scheme.compute_rdf(rdf_grid)
    compute_rdf.assert_called_once_with(scheme.inputs.dimension, rdf_grid)
    scheme.db_handler.insert_results.assert_called_once_with(
        {
            "rdf": scheme.results.rdf,
            "rdf_grid": scheme.results.rdf_grid,
        },
        conflict_mode=DataBaseHandler.ConflictMode.UPDATE,
    )


def test_compute_rdf_without_results(scheme):
    scheme.results = None
    scheme.compute_rdf()
    scheme.db_handler.insert_results.assert_not_called()


def test_input_initialization(mocker):
    dbInfo = mocker.patch("qupled.hf.DatabaseInfo")
    coupling = 1.0
    degeneracy = 2.0
    inputs = hf.Input(coupling, degeneracy)
    assert inputs.coupling == coupling
    assert inputs.degeneracy == degeneracy
    assert inputs.chemical_potential == [-10.0, 10.0]
    assert inputs.cutoff == 10.0
    assert inputs.frequency_cutoff == 10.0
    assert inputs.integral_error == 1.0e-5
    assert inputs.integral_strategy == "full"
    assert inputs.matsubara == 128
    assert inputs.resolution == 0.1
    assert inputs.threads == 1
    assert inputs.theory == "HF"
    assert inputs.database_info == dbInfo.return_value


def test_input_to_native(mocker, inputs):
    native_input = mocker.Mock()
    inputs.to_native(native_input)
    assert native_input.coupling == 1.0
    assert native_input.degeneracy == 2.0


def test_result_initialization(results):
    assert results.idr is None
    assert results.rdf is None
    assert results.rdf_grid is None
    assert results.sdr is None
    assert results.lfc is None
    assert results.ssf is None
    assert results.uint is None
    assert results.wvg is None


def test_result_from_native(mocker, results):
    native_scheme = mocker.Mock()
    native_scheme.idr = np.array([1, 2, 3])
    results.from_native(native_scheme)
    assert np.array_equal(results.idr, np.array([1, 2, 3]))


def test_result_compute_rdf_with_default_grid(mocker, results):
    native_compute_rdf = mocker.patch("qupled.native.compute_rdf")
    mock_dimension = mocker.patch("qupled.native.Dimension")
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.ssf = np.array([4.0, 5.0, 6.0])
    native_compute_rdf.return_value = np.array([7.0, 8.0, 9.0])
    results.compute_rdf("D3")
    assert np.allclose(results.rdf_grid, np.arange(0.0, 10.0, 0.01))
    native_compute_rdf.assert_called_once_with(
        results.rdf_grid, results.wvg, results.ssf, mock_dimension.D3
    )
    assert np.allclose(results.rdf, np.array([7.0, 8.0, 9.0]))


def test_result_compute_rdf_with_custom_grid(mocker, results):
    native_compute_rdf = mocker.patch("qupled.native.compute_rdf")
    mock_dimension = mocker.patch("qupled.native.Dimension")
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.ssf = np.array([4.0, 5.0, 6.0])
    custom_grid = np.array([0.5, 1.5, 2.5])
    native_compute_rdf.return_value = np.array([10.0, 11.0, 12.0])
    results.compute_rdf("D2", custom_grid)
    assert np.allclose(results.rdf_grid, custom_grid)
    native_compute_rdf.assert_called_once_with(
        custom_grid, results.wvg, results.ssf, mock_dimension.D2
    )
    assert np.allclose(results.rdf, np.array([10.0, 11.0, 12.0]))


def test_result_compute_rdf_no_wvg_or_ssf(mocker, results):
    # Should not call native.compute_rdf if wvg or ssf is None
    native_compute_rdf = mocker.patch("qupled.native.compute_rdf")
    results.wvg = None
    results.ssf = np.array([1, 2, 3])
    results.compute_rdf("D3")
    native_compute_rdf.assert_not_called()
    results.wvg = np.array([1, 2, 3])
    results.ssf = None
    results.compute_rdf("D3")
    native_compute_rdf.assert_not_called()


def test_database_info_initialization():
    db_info = hf.DatabaseInfo()
    assert db_info.name == hf.database.DataBaseHandler.DEFAULT_DATABASE_NAME
    assert db_info.run_id is None
    assert db_info.run_table_name == hf.database.DataBaseHandler.RUN_TABLE_NAME


def test_database_info_to_native(mocker):
    native_db_info = mocker.patch("qupled.native.DatabaseInfo")
    db_info = hf.DatabaseInfo()
    db_info.name = "test_db"
    db_info.run_id = 123
    db_info.run_table_name = "test_table"
    native_instance = db_info.to_native()
    assert native_instance == native_db_info.return_value
    assert native_instance.name == "test_db"
    assert native_instance.run_id == 123
    assert native_instance.run_table_name == "test_table"
