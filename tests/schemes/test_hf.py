import pytest
import numpy as np
from unittest.mock import PropertyMock

from qupled import native
from qupled.database.base_tables import ConflictMode, RunStatus
from qupled.database.database_handler import DataBaseHandler
from qupled.schemes import hf
from qupled.util.dimension import Dimension


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


@pytest.mark.unit
def test_native_to_run_status_mapping():
    assert hf.Solver.NATIVE_TO_RUN_STATUS[0] == RunStatus.SUCCESS
    assert hf.Solver.NATIVE_TO_RUN_STATUS[1] == RunStatus.FAILED
    assert len(hf.Solver.NATIVE_TO_RUN_STATUS) == 2


@pytest.mark.unit
def test_hf_initialization():
    scheme = hf.Solver()
    assert scheme.inputs is None
    assert isinstance(scheme.results, hf.Result)
    assert isinstance(scheme.db_handler, DataBaseHandler)
    assert scheme.native_scheme_cls == native.HF
    assert scheme.native_inputs_cls == native.Input
    assert scheme.native_scheme_status is None


@pytest.mark.unit
def test_run_id(scheme):
    run_id = "run_id"
    scheme._db_tables.run_id = run_id
    assert scheme.run_id == run_id


@pytest.mark.unit
def test_compute(scheme, inputs, mocker):
    add_run_to_database = mocker.patch.object(scheme, "_add_run_to_database")
    compute_native = mocker.patch.object(scheme, "_compute_native")
    save = mocker.patch.object(scheme, "_save")
    scheme.compute(inputs)
    assert scheme.inputs is not None
    add_run_to_database.assert_called_once()
    compute_native.assert_called_once()
    save.assert_called_once()


@pytest.mark.unit
def test_db_tables(scheme):
    assert scheme._db_tables is scheme.db_handler.scheme_tables


@pytest.mark.unit
def test_add_run_to_database(scheme, mocker):
    mocker.patch.object(hf.Solver, "run_id", new_callable=PropertyMock).return_value = (
        "mocked-run-id"
    )
    scheme.inputs = mocker.Mock()
    scheme._add_run_to_database()
    scheme._db_tables.insert_run.assert_called_once_with(scheme.inputs)
    assert scheme.inputs.database_info.run_id == scheme.run_id


@pytest.mark.unit
def test_compute_native_with_mpi(scheme, mocker):
    mocker.patch("qupled.native.uses_mpi", True)
    scheme._compute_native_mpi = mocker.Mock()
    scheme._compute_native_serial = mocker.Mock()
    scheme._compute_native()
    scheme._compute_native_mpi.assert_called_once()
    scheme._compute_native_serial.assert_not_called()


@pytest.mark.unit
def test_compute_native_serial(scheme, mocker):
    mocker.patch("qupled.native.uses_mpi", False)
    scheme._compute_native_mpi = mocker.Mock()
    scheme._compute_native_serial = mocker.Mock()
    scheme._compute_native()
    scheme._compute_native_serial.assert_called_once()
    scheme._compute_native_mpi.assert_not_called()


@pytest.mark.unit
def test_compute_native_serial(scheme, inputs, mocker):
    native_input = mocker.Mock()
    native_inputs_cls = mocker.patch.object(
        scheme, "native_inputs_cls", return_value=native_input
    )
    to_native = mocker.patch("qupled.schemes.hf.Input.to_native")
    native_scheme = mocker.Mock()
    native_scheme_cls = mocker.patch.object(
        scheme, "native_scheme_cls", return_value=native_scheme
    )
    from_native = mocker.patch("qupled.schemes.hf.Result.from_native")
    native_scheme.compute.return_value = "mocked-status"
    scheme.inputs = inputs
    scheme._compute_native()
    native_inputs_cls.assert_called_once()
    to_native.assert_called_once_with(native_input)
    native_scheme_cls.assert_called_once_with(native_input)
    native_scheme.compute.assert_called_once()
    assert scheme.native_scheme_status == "mocked-status"
    from_native.assert_called_once_with(native_scheme)


@pytest.mark.unit
def test_compute_native_mpi_calls_all_mpi_functions(scheme, mocker, inputs):
    write_inputs = mocker.patch("qupled.util.mpi.write_inputs")
    launch_mpi_execution = mocker.patch("qupled.util.mpi.launch_mpi_execution")
    read_status = mocker.patch(
        "qupled.util.mpi.read_status", return_value="mocked-status"
    )
    read_results = mocker.patch(
        "qupled.util.mpi.read_results", return_value="mocked-results"
    )
    clean_files = mocker.patch("qupled.util.mpi.clean_files")
    scheme.inputs = inputs
    scheme.results = None
    scheme._compute_native_mpi()
    write_inputs.assert_called_once_with(inputs)
    launch_mpi_execution.assert_called_once_with(scheme.__module__, inputs.processes)
    read_status.assert_called_once()
    read_results.assert_called_once_with(type(None))
    clean_files.assert_called_once()
    assert scheme.native_scheme_status == "mocked-status"
    assert scheme.results == "mocked-results"


@pytest.mark.unit
def test_compute_native_mpi_with_existing_results_type(scheme, mocker, inputs):
    mocker.patch("qupled.util.mpi.write_inputs")
    mocker.patch("qupled.util.mpi.launch_mpi_execution")
    mocker.patch("qupled.util.mpi.read_status", return_value="status")
    read_results = mocker.patch("qupled.util.mpi.read_results", return_value="results")
    mocker.patch("qupled.util.mpi.clean_files")
    scheme.inputs = inputs
    scheme.results = hf.Result()
    scheme._compute_native_mpi()
    read_results.assert_called_once_with(hf.Result)
    assert scheme.results == "results"
    assert scheme.native_scheme_status == "status"


@pytest.mark.unit
def test_run_mpi_worker(mocker):
    mock_inputs = mocker.Mock()
    mock_native_inputs = mocker.Mock()
    mock_scheme = mocker.Mock()
    mock_status = "mocked-status"
    mock_InputCls = mocker.Mock()
    mock_ResultCls = mocker.Mock()
    read_inputs = mocker.patch("qupled.util.mpi.read_inputs", return_value=mock_inputs)
    native_inputs_cls = mocker.patch.object(
        hf.Solver, "native_inputs_cls", return_value=mock_native_inputs
    )
    to_native = mocker.patch.object(mock_inputs, "to_native")
    native_scheme_cls = mocker.patch.object(
        hf.Solver, "native_scheme_cls", return_value=mock_scheme
    )
    mock_scheme.compute.return_value = mock_status
    write_results = mocker.patch("qupled.util.mpi.write_results")
    write_status = mocker.patch("qupled.util.mpi.write_status")
    hf.Solver.run_mpi_worker(mock_InputCls, mock_ResultCls)
    read_inputs.assert_called_once_with(mock_InputCls)
    native_inputs_cls.assert_called_once()
    to_native.assert_called_once_with(mock_native_inputs)
    native_scheme_cls.assert_called_once_with(mock_native_inputs)
    mock_scheme.compute.assert_called_once()
    write_results.assert_called_once_with(mock_scheme, mock_ResultCls)
    write_status.assert_called_once_with(mock_scheme, mock_status)


@pytest.mark.unit
def test_save(scheme, results, mocker):
    scheme.results = results
    scheme.native_scheme_status = mocker.Mock()
    scheme._save()
    scheme._db_tables.update_run_status.assert_called_once_with(
        hf.Solver.NATIVE_TO_RUN_STATUS.get(
            scheme.native_scheme_status, RunStatus.FAILED
        )
    )
    scheme._db_tables.insert_results.assert_called_once_with(scheme.results.__dict__)


@pytest.mark.unit
def test_compute_rdf_with_default_grid(scheme, inputs, results, mocker):
    compute_rdf = mocker.patch("qupled.schemes.hf.Result.compute_rdf")
    scheme.results = results
    scheme.inputs = inputs
    scheme.compute_rdf()
    compute_rdf.assert_called_once_with(scheme.inputs.dimension, None)
    scheme._db_tables.insert_results.assert_called_once_with(
        {
            "rdf": scheme.results.rdf,
            "rdf_grid": scheme.results.rdf_grid,
        },
        conflict_mode=ConflictMode.UPDATE,
    )


@pytest.mark.unit
def test_compute_rdf_with_custom_grid(scheme, inputs, results, mocker):
    compute_rdf = mocker.patch("qupled.schemes.hf.Result.compute_rdf")
    scheme.results = results
    scheme.inputs = inputs
    rdf_grid = np.array([1, 2, 3])
    scheme.compute_rdf(rdf_grid)
    compute_rdf.assert_called_once_with(scheme.inputs.dimension, rdf_grid)
    scheme._db_tables.insert_results.assert_called_once_with(
        {
            "rdf": scheme.results.rdf,
            "rdf_grid": scheme.results.rdf_grid,
        },
        conflict_mode=ConflictMode.UPDATE,
    )


@pytest.mark.unit
def test_compute_rdf_without_results(scheme):
    scheme.results = None
    scheme.compute_rdf()
    scheme.db_handler.insert_results.assert_not_called()


@pytest.mark.unit
def test_input_initialization():
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
    assert inputs.processes == 1
    assert inputs.theory == "HF"
    assert inputs.database_info is None
    assert inputs.dimension == Dimension._3D


@pytest.mark.unit
def test_input_to_native(mocker, inputs):
    native_input = mocker.Mock()
    inputs.to_native(native_input)
    assert native_input.coupling == 1.0
    assert native_input.degeneracy == 2.0


@pytest.mark.unit
def test_result_initialization(results):
    assert results.chemical_potential is None
    assert results.idr is None
    assert results.itcf is None
    assert results.rdf is None
    assert results.rdf_grid is None
    assert results.sdr is None
    assert results.lfc is None
    assert results.ssf is None
    assert results.tau is None
    assert results.uint is None
    assert results.wvg is None


@pytest.mark.unit
def test_result_from_native(mocker, results):
    native_scheme = mocker.Mock()
    native_scheme.idr = np.array([1, 2, 3])
    results.from_native(native_scheme)
    assert np.array_equal(results.idr, np.array([1, 2, 3]))


@pytest.mark.unit
def test_result_compute_rdf_with_default_grid(mocker, results):
    native_compute_rdf = mocker.patch("qupled.native.compute_rdf")
    mock_native_dimension = mocker.patch("qupled.native.Dimension")
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.ssf = np.array([4.0, 5.0, 6.0])
    native_compute_rdf.return_value = np.array([7.0, 8.0, 9.0])
    results.compute_rdf(Dimension._3D)
    assert np.allclose(results.rdf_grid, np.arange(0.0, 10.0, 0.01))
    native_compute_rdf.assert_called_once_with(
        results.rdf_grid, results.wvg, results.ssf, mock_native_dimension.D3
    )
    assert np.allclose(results.rdf, np.array([7.0, 8.0, 9.0]))


@pytest.mark.unit
def test_result_compute_rdf_with_custom_grid(mocker, results):
    native_compute_rdf = mocker.patch("qupled.native.compute_rdf")
    mock_dimension = mocker.patch("qupled.native.Dimension")
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.ssf = np.array([4.0, 5.0, 6.0])
    custom_grid = np.array([0.5, 1.5, 2.5])
    native_compute_rdf.return_value = np.array([10.0, 11.0, 12.0])
    results.compute_rdf(Dimension._2D, custom_grid)
    assert np.allclose(results.rdf_grid, custom_grid)
    native_compute_rdf.assert_called_once_with(
        custom_grid, results.wvg, results.ssf, mock_dimension.D2
    )
    assert np.allclose(results.rdf, np.array([10.0, 11.0, 12.0]))


@pytest.mark.unit
def test_result_compute_rdf_missing_wvg_raises(results):
    results.wvg = None
    results.ssf = np.array([1, 2, 3])
    with pytest.raises(ValueError):
        results.compute_rdf(Dimension._3D)


@pytest.mark.unit
def test_result_compute_rdf_missing_ssf_raises(results):
    results.wvg = np.array([1, 2, 3])
    results.ssf = None
    with pytest.raises(ValueError):
        results.compute_rdf(Dimension._3D)


@pytest.mark.unit
def test_compute_itcf_with_default_grid(scheme, inputs, results, mocker):
    compute_itcf = mocker.patch("qupled.schemes.hf.Result.compute_itcf")
    scheme.results = results
    scheme.inputs = inputs
    scheme.compute_itcf()
    compute_itcf.assert_called_once_with(scheme.inputs, None)
    scheme._db_tables.insert_results.assert_called_once_with(
        {
            "itcf": scheme.results.itcf,
            "tau": scheme.results.tau,
        },
        conflict_mode=ConflictMode.UPDATE,
    )


@pytest.mark.unit
def test_compute_itcf_with_custom_grid(scheme, inputs, results, mocker):
    compute_itcf = mocker.patch("qupled.schemes.hf.Result.compute_itcf")
    scheme.results = results
    scheme.inputs = inputs
    tau = np.array([0.0, 0.1, 0.2])
    scheme.compute_itcf(tau)
    compute_itcf.assert_called_once_with(scheme.inputs, tau)
    scheme._db_tables.insert_results.assert_called_once_with(
        {
            "itcf": scheme.results.itcf,
            "tau": scheme.results.tau,
        },
        conflict_mode=ConflictMode.UPDATE,
    )


@pytest.mark.unit
def test_compute_itcf_without_results(scheme):
    scheme.results = None
    scheme.compute_itcf()
    scheme.db_handler.insert_results.assert_not_called()


@pytest.mark.unit
def test_result_compute_itcf_with_default_grid(mocker, results, inputs):
    # inputs.degeneracy = 2.0, so default tau = np.arange(0.0, 1.1, 0.1) / 2.0
    native_compute_itcf = mocker.patch("qupled.native.compute_itcf_non_interacting")
    native_input = mocker.Mock()
    mocker.patch.object(native, "Input", return_value=native_input)
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.lfc = np.array([4.0, 5.0, 6.0])
    results.chemical_potential = 0.5
    results.idr = np.array([7.0, 8.0, 9.0])
    native_compute_itcf.return_value = np.array([[10.0, 11.0, 12.0]])
    results.compute_itcf(inputs)
    assert np.allclose(results.tau, np.arange(0.0, 1.1, 0.1) / inputs.degeneracy)
    native_compute_itcf.assert_called_once_with(
        native_input,
        results.wvg,
        results.tau,
        results.chemical_potential,
        results.idr,
    )
    assert np.allclose(results.itcf, np.array([[10.0, 11.0, 12.0]]))


@pytest.mark.unit
def test_result_compute_itcf_with_custom_grid(mocker, results, inputs):
    # inputs.degeneracy = 2.0; all custom_tau values are <= 1/2.0, so none filtered
    native_compute_itcf = mocker.patch("qupled.native.compute_itcf_non_interacting")
    native_input = mocker.Mock()
    mocker.patch.object(native, "Input", return_value=native_input)
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.lfc = np.array([4.0, 5.0, 6.0])
    results.chemical_potential = 0.5
    results.idr = np.array([7.0, 8.0, 9.0])
    custom_tau = np.array([0.0, 0.2, 0.4])
    native_compute_itcf.return_value = np.array([[13.0, 14.0, 15.0]])
    results.compute_itcf(inputs, custom_tau)
    assert np.allclose(results.tau, custom_tau)
    native_compute_itcf.assert_called_once_with(
        native_input,
        results.wvg,
        results.tau,
        results.chemical_potential,
        results.idr,
    )
    assert np.allclose(results.itcf, np.array([[13.0, 14.0, 15.0]]))


@pytest.mark.unit
def test_result_compute_itcf_missing_wvg_raises(results, inputs):
    results.wvg = None
    results.lfc = np.array([1, 2, 3])
    with pytest.raises(ValueError):
        results.compute_itcf(inputs)


@pytest.mark.unit
def test_result_compute_itcf_missing_lfc_raises(results, inputs):
    results.wvg = np.array([1, 2, 3])
    results.lfc = None
    with pytest.raises(ValueError):
        results.compute_itcf(inputs)


@pytest.mark.unit
def test_result_compute_itcf_tau_clipped_to_beta(mocker, results, inputs):
    # tau values above 1/theta should be discarded (theta=2.0, beta=0.5)
    mocker.patch("qupled.native.compute_itcf_non_interacting")
    mocker.patch.object(native, "Input")
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.lfc = np.array([4.0, 5.0, 6.0])
    results.chemical_potential = 0.5
    results.idr = np.array([7.0, 8.0, 9.0])
    tau_with_excess = np.array([0.0, 0.2, 0.4, 0.6, 0.8])
    results.compute_itcf(inputs, tau_with_excess)
    assert np.allclose(results.tau, [0.0, 0.2, 0.4])


@pytest.mark.unit
def test_result_invoke_native_itcf(mocker, results, inputs):
    native_compute_itcf = mocker.patch("qupled.native.compute_itcf_non_interacting")
    native_input = mocker.Mock()
    mocker.patch.object(native, "Input", return_value=native_input)
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.tau = np.array([0.0, 0.1, 0.2])
    results.chemical_potential = 0.5
    results.idr = np.array([7.0, 8.0, 9.0])
    expected = np.array([[10.0, 11.0, 12.0]])
    native_compute_itcf.return_value = expected
    result = results._invoke_native_itcf(inputs)
    native_compute_itcf.assert_called_once_with(
        native_input,
        results.wvg,
        results.tau,
        results.chemical_potential,
        results.idr,
    )
    assert np.allclose(result, expected)


@pytest.mark.unit
def test_database_info_initialization():
    db_info = hf.DatabaseInfo()
    assert db_info.blob_storage is None
    assert db_info.name is None
    assert db_info.run_id is None
    assert db_info.run_table_name == hf.SCHEME_RUN_TABLE_NAME


@pytest.mark.unit
def test_database_info_to_native(mocker):
    native_db_info = mocker.patch("qupled.native.DatabaseInfo")
    db_info = hf.DatabaseInfo()
    db_info.blob_storage = "blob_data"
    db_info.name = "test_db"
    db_info.run_id = 123
    db_info.run_table_name = "test_table"
    native_instance = db_info.to_native()
    assert native_instance == native_db_info.return_value
    assert native_instance.blob_storage == "blob_data"
    assert native_instance.name == "test_db"
    assert native_instance.run_id == 123
    assert native_instance.run_table_name == "test_table"
