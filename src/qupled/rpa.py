# -----------------------------------------------------------------------
# RPA class
# -----------------------------------------------------------------------

from __future__ import annotations
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
        scheme = native.Rpa(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)

    # Input class
    class Input:
        """
        Class used to manage the input for the :obj:`qupled.classic.Rpa` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            self.coupling: float = coupling
            """Coupling parameter."""
            self.degeneracy: float = degeneracy
            """Degeneracy parameter."""
            self.chemical_potential: list[float] = [-10.0, 10.0]
            """Initial guess for the chemical potential. Default = ``[-10, 10]``"""
            self.matsubara: int = 128
            """Number of Matsubara frequencies. Default = ``128``"""
            self.resolution: float = 0.1
            """Resolution of the wave-vector grid. Default =  ``0.1``"""
            self.cutoff: float = 10.0
            """Cutoff for the wave-vector grid. Default =  ``10.0``"""
            self.frequency_cutoff: float = 10.0
            """Cutoff for the frequency (applies only in the ground state). Default =  ``10.0``"""
            self.integral_error: float = 1.0e-5
            """Accuracy (relative error) in the computation of integrals. Default = ``1.0e-5``"""
            self.integral_strategy: str = "full"
            """
            Scheme used to solve two-dimensional integrals
            allowed options include:

            - full: the inner integral is evaluated at arbitrary points
              selected automatically by the quadrature rule

            - segregated: the inner integral is evaluated on a fixed
              grid that depends on the integrand that is being processed

            Segregated is usually faster than full but it could become
            less accurate if the fixed points are not chosen correctly. Default =  ``'full'``
            """
            self.threads: int = 1
            """Number of OMP threads for parallel calculations. Default =  ``1``"""
            self.theory: str = "RPA"

        def to_native(self) -> native.RpaInput:
            native_input = native.RpaInput()
            for attr, value in self.__dict__.items():
                setattr(native_input, attr, value)
            return native_input
