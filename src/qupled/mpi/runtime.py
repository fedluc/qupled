from __future__ import annotations

import shutil
import subprocess
import sys

from qupled import native
from qupled.mpi.worker_files import WorkerFiles

# MPI command
MPI_COMMAND = "mpiexec"
_WORKER_MODULE = f"{__package__}.worker"


def _launch_mpi_execution(solver_cls, nproc, worker_files: WorkerFiles):
    """Launch the centralized MPI worker for a solver class.

    Args:
        solver_cls: Solver class to execute in the MPI worker.
        nproc: Number of processes to use for MPI execution.
        worker_files: Files used to exchange data with the MPI worker.

    Raises:
        subprocess.CalledProcessError: If the subprocess execution fails.
    """
    solver_reference = f"{solver_cls.__module__}:{solver_cls.__qualname__}"
    worker_command = [
        sys.executable,
        "-m",
        _WORKER_MODULE,
        "--solver",
        solver_reference,
        "--worker-directory",
        str(worker_files.directory),
    ]
    call_mpi = shutil.which(MPI_COMMAND) is not None and native.uses_mpi
    if call_mpi:
        subprocess.run(
            [MPI_COMMAND, "-n", str(nproc), *worker_command],
            check=True,
        )
    else:
        print("WARNING: Could not call MPI, defaulting to serial execution.")
        subprocess.run(worker_command, check=True)


def run_solver(solver_cls, inputs, nproc, result_cls):
    """Run a solver through the centralized MPI worker.

    Args:
        solver_cls: Solver class to execute in the MPI worker.
        inputs: Serializable input object for the solver run.
        nproc: Number of processes to use for MPI execution.
        result_cls: Serializable result class used to reconstruct results.

    Returns:
        Tuple containing the native solver status and reconstructed results.
    """
    worker_files = WorkerFiles()
    try:
        worker_files.write_inputs(inputs)
        _launch_mpi_execution(solver_cls, nproc, worker_files)
        return worker_files.read_status(), worker_files.read_results(result_cls)
    finally:
        worker_files.cleanup()
