# -----------------------------------------------------------------------
# RPA class
# -----------------------------------------------------------------------

from __future__ import annotations
import numpy as np
from . import native
from . import base


class Rpa(base.ClassicScheme):

    # Compute
    def compute(self, inputs: Rpa.Input):
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        super().compute(inputs, native.Rpa, native.RpaInput(), base.Result())

    # Input class
    class Input(base.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.Rpa` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.theory: str = "RPA"
