# -----------------------------------------------------------------------
# VSStls class
# -----------------------------------------------------------------------

from __future__ import annotations
import pandas as pd
import numpy as np
from . import native
from . import util
from . import stls
from . import base as base


class VSStls(base.IterativeScheme):

    # Compute
    def compute(self, inputs: VSStls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        super().compute(inputs, native.VSStls, self.Result)

    # Set the free energy integrand from a dataframe produced in output
    @staticmethod
    def get_free_energy_integrand(
        run_id: str, database_name: str | None = None
    ) -> native.FreeEnergyIntegrand:
        """Constructs the free energy integrand by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the free energy integrand.
        """
        names = [
            util.HDF.ResultNames.FXC_GRID.value,
            util.HDF.ResultNames.FXC_INT.value,
            util.HDF.ResultNames.ALPHA.value,
        ]
        data = util.HDF.read_results(run_id, database_name, names)
        fxci = native.FreeEnergyIntegrand()
        fxci.grid = data[names[0]]
        fxci.integrand = data[names[1]]
        fxci.alpha = data[names[2]]
        return fxci

    # Input class
    class Input(stls.Stls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.VSStls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.alpha: list[float] = [0.5, 1.0]
            """Initial guess for the free parameter. Default = ``[0.5, 1.0]``"""
            self.coupling_resolution: float = 0.1
            """Resolution of the coupling parameter grid. Default = ``0.1``"""
            self.degeneracy_resolution: float = 0.1
            """Resolution of the degeneracy parameter grid. Default = ``0.1``"""
            self.error_alpha: float = 1.0e-3
            """Minimum error for convergence in the free parameter. Default = ``1.0e-3``"""
            self.iterations_alpha: int = 50
            """Maximum number of iterations to determine the free parameter. Default = ``50``"""
            self.free_energy_integrand: native.FreeEnergyIntegrand = (
                native.FreeEnergyIntegrand()
            )
            """Pre-computed free energy integrand."""
            self.threads: int = 9
            """Number of threads. Default = ``9``"""
            # Undocumented default values
            self.theory: str = "VSSTLS"

        def to_native(self) -> native.VSStlsInput:
            return base.Input.to_native(self, native.VSStlsInput())

    # Results class
    class Result(base.Result):
        """
        Class used to store the results for the :obj:`qupled.classic.VSStls` class.
        """

        def __init__(self, scheme):
            super().__init__()
            self.free_energy_grid = None
            """Free energy grid"""
            self.free_energy_integrand = None
            """Free energy integrand"""
            self.alpha = None
            """Free parameter"""
            super().from_native(scheme)
