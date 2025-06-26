from __future__ import annotations

from . import hf
from . import native
from . import serialize


class Solver(hf.Solver):
    """
    Class used to solve the RPA scheme.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.Rpa

    def __init__(self):
        super().__init__()
        self.results: hf.Result = hf.Result()


@serialize.serializable_dataclass
class Input(hf.Input):
    """
    Class used to manage the input for the :obj:`qupled.rpa.Rpa` class.
    """

    theory: str = "RPA"


if __name__ == "__main__":
    from .mpi_worker import run_mpi_worker

    run_mpi_worker(Solver, Input, hf.Result)
