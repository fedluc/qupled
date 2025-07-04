from __future__ import annotations

from . import hf2d
from . import native
from . import serialize


class Solver(hf2d.Solver):
    """
    Class used to solve the RPA scheme.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.Rpa

    def __init__(self):
        super().__init__()
        self.results: hf2d.Result = hf2d.Result()


@serialize.serializable_dataclass
class Input(hf2d.Input):
    """
    Class used to manage the input for the :obj:`qupled.rpa2d.Rpa2D` class.
    """

    theory: str = "RPA2D"


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, hf2d.Result)