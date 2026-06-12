from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

from qupled import native

# MPI command
MPI_COMMAND = "mpiexec"

# Temporary files used for MPI executions
INPUT_FILE = Path("input.json")
RESULT_FILE = Path("results.json")
STATUS_FILE = Path("status.json")
_WORKER_MODULE = f"{__package__}.worker"


def launch_mpi_execution(solver_cls, nproc):
    """Launch the centralized MPI worker for a solver class.

    Args:
        solver_cls: Solver class to execute in the MPI worker.
        nproc: Number of processes to use for MPI execution.

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
    ]
    call_mpi = shutil.which(MPI_COMMAND) is not None and native.uses_mpi
    if call_mpi:
        subprocess.run([MPI_COMMAND, "-n", str(nproc), *worker_command], check=True)
    else:
        print("WARNING: Could not call MPI, defaulting to serial execution.")
        subprocess.run(worker_command, check=True)


def write_inputs(inputs):
    """Write input data to the MPI input file.

    Args:
        inputs: Serializable input object for the solver run.
    """
    with INPUT_FILE.open("w") as f:
        json.dump(inputs.to_dict(), f)


def read_inputs(input_cls):
    """Read input data from the MPI input file.

    Args:
        input_cls: Serializable input class used to reconstruct the inputs.

    Returns:
        Reconstructed input object.
    """
    with INPUT_FILE.open() as f:
        input_dict = json.load(f)
    return input_cls.from_dict(input_dict)


def write_results(scheme, result_cls):
    """Write solver results from the root MPI rank.

    Args:
        scheme: Native scheme instance.
        result_cls: Serializable result class used to convert native results.
    """
    if scheme.is_root:
        results = result_cls()
        results.from_native(scheme)
        with RESULT_FILE.open("w") as f:
            json.dump(results.to_dict(), f)


def read_results(result_cls):
    """Read solver results from the MPI result file.

    Args:
        result_cls: Serializable result class used to reconstruct the results.

    Returns:
        Reconstructed result object.
    """
    with RESULT_FILE.open() as f:
        result_dict = json.load(f)
    return result_cls.from_dict(result_dict)


def write_status(scheme, status):
    """Write solver status from the root MPI rank.

    Args:
        scheme: Native scheme instance.
        status: Solver status to serialize.
    """
    if scheme.is_root:
        with STATUS_FILE.open("w") as f:
            json.dump(status, f)


def read_status():
    """Read solver status from the MPI status file.

    Returns:
        Deserialized solver status.
    """
    with STATUS_FILE.open() as f:
        status = json.load(f)
    return status


def clean_files():
    """Remove temporary MPI protocol files if they exist."""
    for file in [INPUT_FILE, RESULT_FILE, STATUS_FILE]:
        if file.exists():
            file.unlink()
