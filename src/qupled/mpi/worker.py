from __future__ import annotations

import argparse
import importlib
from collections.abc import Sequence
from pathlib import Path

from qupled.mpi.worker_files import WorkerFiles


def _load_object(reference: str):
    """Load an object from a module-qualified reference.

    Args:
        reference: Object reference in ``module:qualname`` format, such as
            ``qupled.schemes.hf:Solver``.

    Returns:
        Object referenced by ``reference``.

    Raises:
        ValueError: If the reference does not use ``module:qualname`` format.
    """
    # Use ':' to make the module/object boundary explicit.
    module_name, separator, qualname = reference.partition(":")
    if not separator or not module_name or not qualname:
        raise ValueError("Object reference must use 'module:qualname' format.")

    obj = importlib.import_module(module_name)
    for attr in qualname.split("."):
        obj = getattr(obj, attr)
    return obj


def _run_solver_worker(solver_cls, input_cls, result_cls, worker_files: WorkerFiles):
    """Run a native solver inside an MPI worker process.

    Args:
        solver_cls: Solver class to execute.
        input_cls: Serializable input class used to reconstruct inputs.
        result_cls: Serializable result class used to convert native results.
        worker_files: Files used to exchange data with the parent process.
    """
    inputs = worker_files.read_inputs(input_cls)
    native_inputs = solver_cls.native_inputs_cls()
    inputs.to_native(native_inputs)
    scheme = solver_cls.native_scheme_cls(native_inputs)
    status = scheme.compute()
    worker_files.write_results(scheme, result_cls)
    worker_files.write_status(scheme, status)


def run_worker(solver_reference: str, worker_files: WorkerFiles):
    """Run the MPI worker for a solver reference.

    Args:
        solver_reference: Solver reference in ``module:qualname`` format.
        worker_files: Worker files for this MPI worker process.

    Raises:
        ValueError: If the solver does not define MPI input/result metadata.
    """
    solver_cls = _load_object(solver_reference)
    input_cls = getattr(solver_cls, "mpi_input_cls", None)
    result_cls = getattr(solver_cls, "mpi_result_cls", None)
    if input_cls is None or result_cls is None:
        raise ValueError(
            f"{solver_reference} must define mpi_input_cls and mpi_result_cls."
        )

    _run_solver_worker(solver_cls, input_cls, result_cls, worker_files)


def main(argv: Sequence[str] | None = None):
    """Run the MPI worker command-line interface.

    Args:
        argv: Optional argument sequence. Defaults to ``sys.argv``.
    """
    parser = argparse.ArgumentParser(description="Run a qupled MPI worker.")
    parser.add_argument(
        "--solver",
        required=True,
        help="Solver class reference in 'module:qualname' format.",
    )
    parser.add_argument(
        "--worker-directory",
        required=True,
        type=Path,
        help="Directory containing private MPI worker input/output files.",
    )
    args = parser.parse_args(argv)
    run_worker(args.solver, WorkerFiles(args.worker_directory))


if __name__ == "__main__":
    main()
