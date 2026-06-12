from __future__ import annotations

from qupled import native
from qupled.schemes import hf
from qupled.util import serialize


class Solver(hf.Solver):
    """
    Class used to solve the ESA scheme.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.ESA

    def __init__(self):
        super().__init__()
        self.results: hf.Result = hf.Result()
        # Undocumented properties
        self.native_scheme_cls = native.ESA


@serialize.serializable_dataclass
class Input(hf.Input):
    """
    Class used to manage the input for the :obj:`qupled.esa.ESA` class.
    """

    theory: str = "ESA"


# Input and result classes used by the centralized MPI worker.
Solver.mpi_input_cls = Input
Solver.mpi_result_cls = hf.Result
