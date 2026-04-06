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

    def _invoke_native_itcf(self, inputs: Input):
        """
        Invoke the native interacting computation of the imaginary-time
        correlation function (ITCF).

        Args:
            inputs: Input parameters for the ITCF computation.

        Returns:
            np.ndarray: The computed interacting ITCF values.
        """
        native_inputs = native.Input()
        inputs.to_native(native_inputs)
        return native.compute_itcf(
            native_inputs,
            self.wvg,
            self.tau,
            self.chemical_potential,
            self.idr,
            self.lfc,
        )


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, Result)
