from __future__ import annotations

import copy
import numpy as np

from . import mpi
from . import native
from . import output
from . import stls


class VSStls(stls.Stls):
    """
    Class used to solve the VSStls scheme.
    """

    def __init__(self):
        super().__init__()
        self.results: Result = Result()
        # Undocumented properties
        self.native_scheme_cls = native.VSStls
        self.native_inputs = native.VSStlsInput()

    def compute(self, inputs: Input):
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self._fill_free_energy_integrand(inputs)
        super().compute(inputs)

    def _fill_free_energy_integrand(self, inputs: Input):
        """
        Fills the free energy integrand by computing missing state points for a given input.

        This method iterates over the missing state points (coupling values) for the given input,
        computes the required data for each coupling value, and updates the input data accordingly.
        After processing, it restores the original coupling value in the input.

        Args:
            inputs: Input parameters.

        Behavior:
            - Identifies the missing state points (coupling values) that need to be computed.
            - For each missing coupling value:
                - Prints a message indicating the current computation.
                - Updates the input coupling value to the current coupling.
                - Performs the computation using the `compute` method.
                - Updates the input data with the results of the computation.
            - Restores the original target coupling value in the input after processing.
        """
        target_coupling = inputs.coupling
        for coupling in self._get_missing_state_points(inputs):
            print("---------------------------------------------------------------")
            print(f"Subcall: solving {inputs.theory} scheme for rs = {coupling:.5f}")
            inputs.coupling = coupling
            self.compute(inputs)
            self._update_input_data(inputs)
        inputs.coupling = target_coupling

    @staticmethod
    def _get_missing_state_points(inputs: Input) -> np.ndarray:
        """
        Calculate the missing state points in a grid based on the expected and actual values.

        This function determines the points in the expected grid that are not present
        in the actual grid of the free energy integrand. The precision of the comparison
        is determined by the resolution of the coupling parameter.

        Args:
            inputs (Input): The input parameters.

        Returns:
            np.ndarray: An array of missing state points in the grid. If the actual grid
            is `None`, the function returns the entire expected grid.
        """
        rs = inputs.coupling
        drs = inputs.coupling_resolution
        expected_grid = np.arange(drs, rs - 0.1 * drs, 3 * drs)
        actual_grid = inputs.free_energy_integrand.grid
        precision = int(np.abs(np.log10(drs)))
        return (
            np.setdiff1d(
                np.round(expected_grid, precision), np.round(actual_grid, precision)
            )
            if actual_grid is not None
            else expected_grid
        )

    def _update_input_data(self, inputs: Input):
        """
        Updates the input data with a free energy integrand object.

        This method creates a `FreeEnergyIntegrand` instance using the results
        stored in the current object and assigns it to the `free_energy_integrand`
        attribute of the provided `inputs` object.

        Args:
            inputs: Input parameters.
        """
        free_energy_integrand = FreeEnergyIntegrand(
            self.results.free_energy_grid,
            self.results.free_energy_integrand,
            self.results.alpha,
        )
        inputs.free_energy_integrand = free_energy_integrand

    # Get the free energy integrand from database
    @staticmethod
    def get_free_energy_integrand(
        run_id: int, database_name: str | None = None
    ) -> FreeEnergyIntegrand:
        """
        Retrieve the free energy integrand for a given run ID from the database.

        Args:
            run_id: The unique identifier for the run whose data is to be retrieved.
            database_name: The name of the database to query.
                If None, the default database will be used.

        Returns:
            native.FreeEnergyIntegrand: An object containing the free energy grid,
            integrand, and alpha values retrieved from the database.
        """
        names = ["free_energy_grid", "free_energy_integrand", "alpha"]
        data = output.DataBase.read_results(run_id, database_name, names)
        return FreeEnergyIntegrand(data[names[0]], data[names[1]], data[names[2]])


# Input class
class Input(stls.Input):
    """
    Class used to manage the input for the :obj:`qupled.vsstls.VSStls` class.
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
        self.free_energy_integrand: FreeEnergyIntegrand = FreeEnergyIntegrand()
        """Pre-computed free energy integrand."""
        self.threads: int = 9
        """Number of threads. Default = ``9``"""
        # Undocumented default values
        self.theory: str = "VSSTLS"


class Result(stls.Result):
    """
    Class used to store the results for the :obj:`qupled.vsstls.VSStls` class.
    """

    def __init__(self):
        super().__init__()
        self.free_energy_grid = None
        """Free energy grid"""
        self.free_energy_integrand = None
        """Free energy integrand"""
        self.alpha = None
        """Free parameter"""


class FreeEnergyIntegrand:

    def __init__(
        self,
        grid: np.ndarray | None = None,
        integrand: np.ndarray | None = None,
        alpha: np.ndarray | None = None,
    ):
        self.grid = grid
        """ Coupling parameter grid. Default = ``None``"""
        self.integrand = integrand
        """ Free energy integrand. Default = ``None``"""
        self.alpha = alpha
        """ Free parameter. Default = ``None``"""

    def to_native(self) -> native.FreeEnergyIntegrand:
        """
        Converts the current object to a native `FreeEnergyIntegrand` instance.

        This method creates an instance of `native.FreeEnergyIntegrand` and maps
        the attributes of the current object to the corresponding attributes of
        the native instance. If an attribute's value is `None`, it is replaced
        with an empty NumPy array.

        Returns:
            native.FreeEnergyIntegrand: A new instance of `FreeEnergyIntegrand`
            with attributes copied from the current object.
        """
        native_guess = native.FreeEnergyIntegrand()
        for attr, value in self.__dict__.items():
            if value is not None:
                setattr(native_guess, attr, value)
        return native_guess
