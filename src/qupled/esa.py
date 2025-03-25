# -----------------------------------------------------------------------
# ESA class
# -----------------------------------------------------------------------

from __future__ import annotations
from . import native
from . import rpa
from . import base


class ESA(base.ClassicScheme):

    # Compute
    def compute(self, inputs: ESA.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        super().compute(inputs, native.ESA, native.RpaInput(), base.Result())

    # Input class
    class Input(rpa.Rpa.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.ESA` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            # Undocumented default values
            self.theory = "ESA"
