from __future__ import annotations

from dataclasses import dataclass

from . import hf
from . import native


class Rpa(hf.HF):
    """
    Class used to solve the RPA scheme.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.Rpa

    def __init__(self):
        super().__init__()
        self.results: hf.Result = hf.Result()


@dataclass
class Input(hf.Input):
    """
    Class used to manage the input for the :obj:`qupled.rpa.Rpa` class.
    """

    theory: str = "RPA"


if __name__ == "__main__":
    from .mpi_worker import run_mpi_worker

    run_mpi_worker(Input, hf.Result, Rpa.native_inputs_cls, Rpa.native_scheme_cls)
