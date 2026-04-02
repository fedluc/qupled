from __future__ import annotations

import numpy as np

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
        self.results: Result = Result()


@serialize.serializable_dataclass
class Input(hf.Input):
    """
    Class used to manage the input for the :obj:`qupled.rpa.Rpa` class.
    """

    theory: str = "RPA"


@serialize.serializable_dataclass
class Result(hf.Result):
    """
    Class used to store the results for the :obj:`qupled.rpa.Solver` class.
    """

    def compute_itcf(self, inputs: Input, tau: np.ndarray | None = None):
        """
        Compute the imaginary-time correlation function (ITCF) for the system.

        Args:
            tau (np.ndarray | None, optional): A 1D array specifying the imaginary-time points
                at which the ITCF is computed. If None, a default grid ranging from 0.0
                to 10.0 with a step size of 0.01 is used.

        Returns:
            None: The computed ITCF is stored in the `self.itcf` attribute.
        """
        if self.wvg is not None and self.lfc is not None:
            self.tau = tau if tau is not None else np.arange(0.0, 0.6, 0.1)
            native_inputs = native.Input()
            inputs.to_native(native_inputs)
            self.itcf = native.compute_itcf(
                native_inputs,
                self.wvg,
                self.tau,
                self.chemical_potential,
                self.idr,
                self.lfc,
            )


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, Result)
