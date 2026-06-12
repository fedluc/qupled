import pytest

from qupled.mpi import worker

pytestmark = pytest.mark.unit


class WorkerInput:
    pass


class WorkerResult:
    pass


class WorkerSolver:
    mpi_input_cls = WorkerInput
    mpi_result_cls = WorkerResult


class MissingMetadataSolver:
    pass


def test_load_object():
    reference = f"{WorkerSolver.__module__}:{WorkerSolver.__qualname__}"
    assert worker._load_object(reference) is WorkerSolver


def test_load_object_rejects_invalid_reference():
    with pytest.raises(ValueError, match="module:qualname"):
        worker._load_object("NotAnImportableObject")


def test_run_worker_uses_solver_mpi_metadata(mocker):
    worker_files = mocker.Mock()
    run_solver_worker = mocker.patch("qupled.mpi.worker._run_solver_worker")
    worker.run_worker(
        f"{WorkerSolver.__module__}:{WorkerSolver.__qualname__}", worker_files
    )
    run_solver_worker.assert_called_once_with(
        WorkerSolver, WorkerInput, WorkerResult, worker_files
    )


def test_run_worker_rejects_missing_mpi_metadata(mocker):
    with pytest.raises(ValueError, match="mpi_input_cls and mpi_result_cls"):
        worker.run_worker(
            f"{MissingMetadataSolver.__module__}:{MissingMetadataSolver.__qualname__}",
            mocker.Mock(),
        )


def test_main_runs_worker(mocker, tmp_path):
    run_worker = mocker.patch("qupled.mpi.worker.run_worker")
    worker.main(
        [
            "--solver",
            "some.module:Solver",
            "--worker-directory",
            str(tmp_path),
        ]
    )
    run_worker.assert_called_once()
    solver_reference, worker_files = run_worker.call_args.args
    assert solver_reference == "some.module:Solver"
    assert worker_files.directory == tmp_path


def test_run_solver_worker_executes_native_solver(mocker):
    mock_inputs = mocker.Mock()
    mock_native_inputs = mocker.Mock()
    mock_scheme = mocker.Mock()
    mock_status = "mocked-status"
    mock_InputCls = mocker.Mock()
    mock_ResultCls = mocker.Mock()
    worker_files = mocker.Mock()
    worker_files.read_inputs.return_value = mock_inputs
    native_inputs_cls = mocker.patch.object(
        WorkerSolver, "native_inputs_cls", return_value=mock_native_inputs, create=True
    )
    to_native = mocker.patch.object(mock_inputs, "to_native")
    native_scheme_cls = mocker.patch.object(
        WorkerSolver, "native_scheme_cls", return_value=mock_scheme, create=True
    )
    mock_scheme.compute.return_value = mock_status
    worker._run_solver_worker(WorkerSolver, mock_InputCls, mock_ResultCls, worker_files)
    worker_files.read_inputs.assert_called_once_with(mock_InputCls)
    native_inputs_cls.assert_called_once_with()
    to_native.assert_called_once_with(mock_native_inputs)
    native_scheme_cls.assert_called_once_with(mock_native_inputs)
    mock_scheme.compute.assert_called_once_with()
    worker_files.write_results.assert_called_once_with(mock_scheme, mock_ResultCls)
    worker_files.write_status.assert_called_once_with(mock_scheme, mock_status)
