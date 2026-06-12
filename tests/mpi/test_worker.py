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
    def run_mpi_worker(cls, input_cls, result_cls):
        raise AssertionError("run_mpi_worker should be patched in tests")


class MissingMetadataSolver:
    @classmethod
    def run_mpi_worker(cls, input_cls, result_cls):
        raise AssertionError("run_mpi_worker should not be called")


def test_load_object():
    reference = f"{WorkerSolver.__module__}:{WorkerSolver.__qualname__}"
    assert worker._load_object(reference) is WorkerSolver


def test_load_object_rejects_invalid_reference():
    with pytest.raises(ValueError, match="module:qualname"):
        worker._load_object("NotAnImportableObject")


def test_run_solver_worker_uses_solver_mpi_metadata(mocker):
    run_mpi_worker = mocker.patch.object(WorkerSolver, "run_mpi_worker")
    worker.run_solver_worker(f"{WorkerSolver.__module__}:{WorkerSolver.__qualname__}")
    run_mpi_worker.assert_called_once_with(WorkerInput, WorkerResult)


def test_run_solver_worker_rejects_missing_mpi_metadata():
    with pytest.raises(ValueError, match="mpi_input_cls and mpi_result_cls"):
        worker.run_solver_worker(
            f"{MissingMetadataSolver.__module__}:{MissingMetadataSolver.__qualname__}"
        )


def test_main_runs_solver_worker(mocker):
    run_solver_worker = mocker.patch("qupled.mpi.worker.run_solver_worker")
    worker.main(["--solver", "some.module:Solver"])
    run_solver_worker.assert_called_once_with("some.module:Solver")
