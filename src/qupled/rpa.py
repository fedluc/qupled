# -----------------------------------------------------------------------
# RPA class
# -----------------------------------------------------------------------

from __future__ import annotations
import numpy as np
from . import native
from . import util
from . import base


class Rpa(base.ClassicScheme):

    INPUT_TABLE_NAME = "RPA_input"
    RESULT_TABLE_NAME = "RPA_results"

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
        results = self.Results(scheme)
        db_handler = base.DataBaseHandler(
            inputs, results, self.INPUT_TABLE_NAME, self.RESULT_TABLE_NAME
        )
        db_handler.insert()

    # Input class
    class Input:
        """
        Class used to manage the input for the :obj:`qupled.classic.Rpa` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            self.chemical_potential: list[float] = [-10.0, 10.0]
            """Initial guess for the chemical potential. Default = ``[-10, 10]``"""
            self.coupling: float = coupling
            """Coupling parameter."""
            self.cutoff: float = 10.0
            """Cutoff for the wave-vector grid. Default =  ``10.0``"""
            self.degeneracy: float = degeneracy
            """Degeneracy parameter."""
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
            self.matsubara: int = 128
            """Number of Matsubara frequencies. Default = ``128``"""
            self.resolution: float = 0.1
            """Resolution of the wave-vector grid. Default =  ``0.1``"""
            self.threads: int = 1
            """Number of OMP threads for parallel calculations. Default =  ``1``"""
            self.theory: str = "RPA"

        def to_native(self) -> native.RpaInput:
            return base.Input.to_native(self, native.RpaInput())

    # Results class
    class Results:
        """
        Class used to store the results for the :obj:`qupled.classic.Rpa` class.
        """

        def __init__(self, scheme):
            self.idr: np.ndarray = None
            """Ideal density response"""
            self.sdr: np.ndarray = None
            """Static density response"""
            self.slfc: np.ndarray = None
            """Static local field correction"""
            self.ssf: np.ndarray = None
            """Static structure factor"""
            self.wvg: np.ndarray = None
            """Wave-vector grid"""
            base.Result.from_native(self, scheme)
