from __future__ import annotations

import numpy as np

from qupled import native
from qupled.schemes import hf
from qupled.util import serialize


class Solver(hf.Solver):
    """
    Class used to solve the reconstruction-driven STLS scheme.
    """

    native_scheme_cls = native.RecStls
    native_inputs_cls = native.RecStlsInput

    def __init__(self):
        super().__init__()
        self.results: Result = Result()
        self.native_scheme_cls = native.RecStls


@serialize.serializable_dataclass
class Input(hf.Input):
    """
    Input for the :obj:`qupled.recstls.RecStls` class.
    """

    rdf_grid: np.ndarray = None
    rdf: np.ndarray = None
    theory: str = "REC-STLS"


@serialize.serializable_dataclass
class Result(hf.Result):
    """
    Results for the :obj:`qupled.recstls.RecStls` class.
    """

    ssf_input: np.ndarray = None


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, Result)
