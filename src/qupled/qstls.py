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
    @util.MPI.record_time
    @util.MPI.synchronize_ranks
    def compute(self, inputs: Qstls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.Qstls(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)

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

        def to_native(self) -> native.QstlsInput:
            native_input = native.QstlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.to_native())
                else:
                    setattr(native_input, attr, value)
            return native_input
