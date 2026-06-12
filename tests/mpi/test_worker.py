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

    @classmethod
    def run_mpi_worker(cls, input_cls, result_cls, worker_files):
        raise AssertionError("run_mpi_worker should be patched in tests")


class MissingMetadataSolver:
    @classmethod
    def run_mpi_worker(cls, input_cls, result_cls, worker_files):
        raise AssertionError("run_mpi_worker should not be called")


def test_load_object():
    reference = f"{WorkerSolver.__module__}:{WorkerSolver.__qualname__}"
    assert worker._load_object(reference) is WorkerSolver


def test_load_object_rejects_invalid_reference():
    with pytest.raises(ValueError, match="module:qualname"):
        worker._load_object("NotAnImportableObject")


def test_run_solver_worker_uses_solver_mpi_metadata(mocker):
    worker_files = mocker.Mock()
    run_mpi_worker = mocker.patch.object(WorkerSolver, "run_mpi_worker")
    worker.run_solver_worker(
        f"{WorkerSolver.__module__}:{WorkerSolver.__qualname__}", worker_files
    )
    run_mpi_worker.assert_called_once_with(WorkerInput, WorkerResult, worker_files)


def test_run_solver_worker_rejects_missing_mpi_metadata(mocker):
    with pytest.raises(ValueError, match="mpi_input_cls and mpi_result_cls"):
        worker.run_solver_worker(
            f"{MissingMetadataSolver.__module__}:{MissingMetadataSolver.__qualname__}",
            mocker.Mock(),
        )


def test_main_runs_solver_worker(mocker, tmp_path):
    run_solver_worker = mocker.patch("qupled.mpi.worker.run_solver_worker")
    worker.main(
        [
            "--solver",
            "some.module:Solver",
            "--worker-directory",
            str(tmp_path),
        ]
    )
    run_solver_worker.assert_called_once()
    solver_reference, worker_files = run_solver_worker.call_args.args
    assert solver_reference == "some.module:Solver"
    assert worker_files.directory == tmp_path
