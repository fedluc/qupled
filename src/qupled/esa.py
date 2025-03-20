# -----------------------------------------------------------------------
# ESA class
# -----------------------------------------------------------------------

from __future__ import annotations
from . import native
from . import util
from . import rpa
from . import base


class ESA(base.ClassicScheme):

    INPUT_TABLE_NAME = "ESA_input"
    RESULT_TABLE_NAME = "ESA_results"

    # Compute
    @util.MPI.record_time
    @util.MPI.synchronize_ranks
    def compute(self, inputs: ESA.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.ESA(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)
        self._save(
            scheme,
            rpa.Rpa.Results(scheme),
            self.INPUT_TABLE_NAME,
            self.RESULT_TABLE_NAME,
        )

    # Input class
    class Input(rpa.Rpa.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.ESA` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            # Undocumented default values
            self.theory = "ESA"
