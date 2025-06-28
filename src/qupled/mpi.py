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

def mpi_process(func):
    def wrapper(*args, **kwargs):
        native.MPI.init()
        try:
            return func(*args, **kwargs)
        finally:
            native.MPI.finalize()

    return wrapper


def launch_mpi_execution(module, nproc):
    call_mpi = shutil.which(MPI_COMMAND) is not None and native.MPI.is_used
    if call_mpi:
        subprocess.run([MPI_COMMAND, "-n", str(nproc), "python", "-m", module], check=True)
    else:
        print("WARNING: Could not call MPI, defaulting to serial execution.")
        subprocess.run(["python", "-m", module], check=True)


def write_inputs(inputs):
    """
    Writes the provided inputs to a file in JSON format.

    Args:
        inputs: An object with a `to_dict()` method that returns a dictionary representation of the inputs.

    Raises:
        IOError: If there is an error opening or writing to the file.
    """
    with INPUT_FILE.open("w") as f:
        json.dump(inputs.to_dict(), f)

def read_inputs(InputCls):
    """
    Load inputs from a JSON file.

    Args:
        file_path (str): Path to the JSON file containing input data.

    Returns:
        InputCls: An instance of the InputCls class populated with data from the file.
    """
    with INPUT_FILE.open() as f:
        input_dict = json.load(f)
    return InputCls.from_dict(input_dict)


def write_results(scheme, ResultCls):
    """
    Write results to a JSON file.

    Args:
        results: The results object to be written to the file.
        file_path (str): Path to the JSON file where results will be saved.
    """
    if scheme.is_root:
        results = ResultCls()
        results.from_native(scheme)
        with RESULT_FILE.open("w") as f:
            json.dump(results.to_dict(), f)

def read_results(ResultsCls):
    """
    Loads results from a JSON file and returns an instance of the specified ResultsCls.

    Args:
        ResultsCls (type): The class with a `from_dict` method to instantiate from the loaded dictionary.

    Returns:
        An instance of ResultsCls initialized with data loaded from the JSON file.

    Raises:
        FileNotFoundError: If the result file does not exist.
        json.JSONDecodeError: If the file content is not valid JSON.
        AttributeError: If ResultsCls does not have a `from_dict` method.
    """
    with RESULT_FILE.open() as f:
        result_dict = json.load(f)
    return ResultsCls.from_dict(result_dict)