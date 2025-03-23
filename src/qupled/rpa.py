# -----------------------------------------------------------------------
# RPA class
# -----------------------------------------------------------------------

from __future__ import annotations
import numpy as np
from . import native
from . import util
from . import base


class Rpa(base.ClassicScheme):

    # Compute
    @util.MPI.record_time
    @util.MPI.synchronize_ranks
    def compute(self, inputs: Rpa.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self.inputs = inputs
        scheme = native.Rpa(self.inputs.to_native())
        self._compute(scheme)
        self.results = self.Result(scheme)
        self._save(scheme)

    # Input class
    class Input(base.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.Rpa` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.theory: str = "RPA"

        def to_native(self) -> native.RpaInput:
            return super().to_native(native.RpaInput())

    # Results class
    class Result(base.Result):
        """
        Class used to store the results for the :obj:`qupled.classic.Rpa` class.
        """

        def __init__(self, scheme):
            super().__init__()
            super().from_native(scheme)
