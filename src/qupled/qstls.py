# -----------------------------------------------------------------------
# Qstls class
# -----------------------------------------------------------------------

from __future__ import annotations
from . import native
from . import util
from . import base
from . import stls


class Qstls(base.QuantumIterativeScheme):

    # Compute
    def compute(self, inputs: Qstls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        super().compute(inputs, native.Qstls, native.QstlsInput(), self.Result())

    # Input class
    class Input(stls.Stls.Input):
        """
        Class used to manage the input for the :obj:`qupled.quantum.Qstls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.fixed: str = ""
            """ Name of the file storing the fixed component of the auxiliary density 
            response in the QSTLS scheme. """
            self.guess: Qstls.Guess = Qstls.Guess()
            """Initial guess. Default = ``Qstls.Guess()``"""
            # Undocumented default values
            self.theory = "QSTLS"

    # Results class
    class Result(base.Result):
        """
        Class used to store the results for the :obj:`qupled.classic.VSStls` class.
        """

        def __init__(self):
            super().__init__()
            self.adr = None
            """Auxiliary density response"""
