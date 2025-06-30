import json
import subprocess
import shutil

from pathlib import Path
from . import native

# MPI command
MPI_COMMAND = "mpiexec"

# Temporary files used for MPI executions
INPUT_FILE = Path("input.json")
RESULT_FILE = Path("results.json")
STATUS_FILE = Path("status.json")


def launch_mpi_execution(module, nproc):
    """
    Launches the execution of a Python module using MPI if available, otherwise defaults to serial execution.

    Args:
        module (str): The name of the Python module to execute (as used with the '-m' flag).
        nproc (int): The number of processes to use for MPI execution.

    Behavior:
        - Checks if the MPI command is available and if native MPI usage is enabled.
        - If MPI is available, runs the module with the specified number of processes using MPI.
        - If MPI is not available, prints a warning and runs the module in serial mode.

    Raises:
        subprocess.CalledProcessError: If the subprocess execution fails.
    """
    call_mpi = shutil.which(MPI_COMMAND) is not None and native.uses_mpi
    if call_mpi:
        subprocess.run(
            [MPI_COMMAND, "-n", str(nproc), "python", "-m", module], check=True
        )
    else:
        print("WARNING: Could not call MPI, defaulting to serial execution.")
        subprocess.run(["python", "-m", module], check=True)


def write_inputs(inputs):
    """
    Writes the input data to the INPUT_FILE in JSON format.
    """
    with INPUT_FILE.open("w") as f:
        json.dump(inputs.to_dict(), f)


def read_inputs(InputCls):
    """
    Reads input data from a predefined input file and constructs an instance of the specified input class.
    """
    with INPUT_FILE.open() as f:
        input_dict = json.load(f)
    return InputCls.from_dict(input_dict)


def write_results(scheme, ResultCls):
    """
    Writes the results of a computation to a JSON file if the current process is the root.
    """
    if scheme.is_root:
        results = ResultCls()
        results.from_native(scheme)
        with RESULT_FILE.open("w") as f:
            json.dump(results.to_dict(), f)


def read_results(ResultsCls):
    """
    Reads results from a JSON file and returns an instance of the specified ResultsCls.
    """
    with RESULT_FILE.open() as f:
        result_dict = json.load(f)
    return ResultsCls.from_dict(result_dict)


def write_status(scheme, status):
    """
    Writes the status of a computation to a JSON file if the current process is the root.
    """
    if scheme.is_root:
        with STATUS_FILE.open("w") as f:
            json.dump(status, f)


def read_status():
    """
    Reads status from a JSON file and returns an instance of the specified ResultsCls.
    """
    with STATUS_FILE.open() as f:
        status = json.load(f)
    return status


def clean_files():
    """
    Removes the input and result files if they exist.
    """
    for file in [INPUT_FILE, RESULT_FILE, STATUS_FILE]:
        if file.exists():
            file.unlink()
