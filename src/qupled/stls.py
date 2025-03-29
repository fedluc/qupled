# -----------------------------------------------------------------------
# Stls class
# -----------------------------------------------------------------------

from __future__ import annotations
from . import native
from . import base


class Stls(base.IterativeScheme):

    # Compute
    def compute(self, inputs: Stls.Input):
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        super().compute(inputs, native.Stls, native.StlsInput(), self.Result())

    # Input class
    class Input(base.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.Stls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.error: float = 1.0e-5
            """Minimum error for convergence. Default = ``1.0e-5``"""
            self.mixing: float = 1.0
            """Mixing parameter. Default = ``1.0``"""
            self.iterations: int = 1000
            """Maximum number of iterations. Default = ``1000``"""
            self.output_frequency: int = 10
            """Output frequency to write the recovery file. Default = ``10``"""
            self.recovery_file: str = ""
            """Name of the recovery file. Default = ``""``"""
            self.guess: Stls.Guess = Stls.Guess()
            """Initial guess. Default = ``Stls.Guess()``"""
            # Undocumented default values
            self.theory: str = "STLS"

    # Results class
    class Result(base.Result):
        """
        Class used to store the results for the :obj:`qupled.classic.Stls` class.
        """

        def __init__(self):
            super().__init__()
            self.error: float = None
            """Residual error in the solution"""
