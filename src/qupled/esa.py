from __future__ import annotations

from dataclasses import dataclass

from . import hf
from . import native


class ESA(hf.HF):
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


@dataclass
class Input(hf.Input):
    """
    Class used to manage the input for the :obj:`qupled.esa.ESA` class.
    """

    theory: str = "ESA"


if __name__ == "__main__":
    from .mpi_worker import run_mpi_worker

    run_mpi_worker(Input, hf.Result, ESA.native_inputs_cls, ESA.native_scheme_cls)
