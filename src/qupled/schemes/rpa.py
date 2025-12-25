from __future__ import annotations

from qupled import native
from qupled.schemes import hf
from qupled.util import serialize


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
    Solver.run_mpi_worker(Input, hf.Result)
