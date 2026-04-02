from __future__ import annotations

from qupled import native
from qupled.schemes import hf, rpa
from qupled.util import serialize


class Solver(hf.Solver):
    """
    Class used to solve the ESA scheme.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.ESA

    def __init__(self):
        super().__init__()
        self.results: rpa.Result = rpa.Result()
        # Undocumented properties
        self.native_scheme_cls = native.ESA


@serialize.serializable_dataclass
class Input(hf.Input):
    """
    Class used to manage the input for the :obj:`qupled.esa.ESA` class.
    """

    theory: str = "ESA"


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, rpa.Result)
