import json
from pathlib import Path


def run_mpi_worker(SolverCls, InputCls, ResultCls):
    inputs = load_inputs_from_file("input.json", InputCls)
    native_inputs = SolverCls.native_inputs_cls()
    inputs.to_native(native_inputs)
    scheme = SolverCls.native_scheme_cls(native_inputs)
    scheme.compute()
    write_results_to_file("results.json", scheme, ResultCls)


def load_inputs_from_file(file_path: str, InputCls):
    """
    Load inputs from a JSON file.

    Args:
        file_path (str): Path to the JSON file containing input data.

    Returns:
        InputCls: An instance of the InputCls class populated with data from the file.
    """
    with Path(file_path).open() as f:
        input_dict = json.load(f)
    return InputCls.from_dict(input_dict)


def write_results_to_file(file_path: str, scheme, ResultCls):
    """
    Write results to a JSON file.

    Args:
        results: The results object to be written to the file.
        file_path (str): Path to the JSON file where results will be saved.
    """
    if scheme.is_root:
        results = ResultCls()
        results.from_native(scheme)
        with Path(file_path).open("w") as f:
            json.dump(results.to_dict(), f)
