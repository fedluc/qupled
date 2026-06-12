from __future__ import annotations

import argparse
import importlib
from collections.abc import Sequence


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


def run_solver_worker(solver_reference: str):
    """Run the MPI worker for a solver reference.

    Args:
        solver_reference: Solver reference in ``module:qualname`` format.

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

    solver_cls.run_mpi_worker(input_cls, result_cls)


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
    args = parser.parse_args(argv)
    run_solver_worker(args.solver)


if __name__ == "__main__":
    main()
