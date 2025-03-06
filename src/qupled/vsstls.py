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
    @util.MPI.record_time
    @util.MPI.synchronize_ranks
    def compute(self, inputs: VSStls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.VSStls(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)

    # Save results
    @util.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        pd.DataFrame(scheme.free_energy_grid).to_hdf(
            self.hdf_file_name, key=util.HDF.EntryKeys.FXC_GRID.value
        )
        pd.DataFrame(scheme.free_energy_integrand).to_hdf(
            self.hdf_file_name, key=util.HDF.EntryKeys.FXC_INT.value
        )
        pd.DataFrame(scheme.alpha).to_hdf(
            self.hdf_file_name, key=util.HDF.EntryKeys.ALPHA.value
        )

    # Set the free energy integrand from a dataframe produced in output
    @staticmethod
    def get_free_energy_integrand(file_name: str) -> native.FreeEnergyIntegrand:
        """Constructs the free energy integrand by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the free energy integrand.
        """
        fxci = native.FreeEnergyIntegrand()
        hdf_data = util.HDF().read(
            file_name,
            [
                util.HDF.EntryKeys.FXC_GRID.value,
                util.HDF.EntryKeys.FXC_INT.value,
                util.HDF.EntryKeys.ALPHA.value,
            ],
        )
        fxci.grid = hdf_data[util.HDF.EntryKeys.FXC_GRID.value]
        fxci.integrand = np.ascontiguousarray(
            hdf_data[util.HDF.EntryKeys.FXC_INT.value]
        )
        fxci.alpha = hdf_data[util.HDF.EntryKeys.ALPHA.value]
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
            native_input = native.VSStlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.to_native())
                else:
                    setattr(native_input, attr, value)
            return native_input
