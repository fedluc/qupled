from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

from . import hf
from . import native
from . import output
from . import stls


class VSStls(stls.Stls):
    """
    Class used to solve the VSStls scheme.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.VSStls
    native_inputs_cls = native.VSStlsInput

    def __init__(self):
        super().__init__()
        self.results: Result = Result()

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
        missing_state_points = self._get_missing_state_points(inputs)
        do_subcall = len(missing_state_points) > 0
        for coupling in missing_state_points:
            print("---------------------------------------------------------------")
            print(f"Subcall: solving {inputs.theory} scheme for rs = {coupling:.5f}")
            inputs.coupling = coupling
            self.compute(inputs)
            self._update_input_data(inputs)
        if do_subcall:
            print("---------------------------------------------------------------")
            print("Subcalls completed.")
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
            self.results.free_energy_grid, self.results.free_energy_integrand
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
            and integrand values retrieved from the database.
        """
        names = ["free_energy_grid", "free_energy_integrand"]
        data = output.DataBase.read_results(run_id, database_name, names)
        return FreeEnergyIntegrand(data[names[0]], data[names[1]])


@dataclass
class Input(stls.Input):
    """
    Class used to manage the input for the :obj:`qupled.vsstls.VSStls` class.
    """

    alpha: list[float] = field(default_factory=lambda: [0.5, 1.0])
    """Initial guess for the free parameter. Default = ``[0.5, 1.0]``"""
    coupling_resolution: float = 0.1
    """Resolution of the coupling parameter grid. Default = ``0.1``"""
    degeneracy_resolution: float = 0.1
    """Resolution of the degeneracy parameter grid. Default = ``0.1``"""
    error_alpha: float = 1.0e-3
    """Minimum error for convergence in the free parameter. Default = ``1.0e-3``"""
    iterations_alpha: int = 50
    """Maximum number of iterations to determine the free parameter. Default = ``50``"""
    free_energy_integrand: FreeEnergyIntegrand = field(
        default_factory=lambda: FreeEnergyIntegrand()
    )
    """Pre-computed free energy integrand."""
    threads: int = 9
    """Number of threads. Default = ``9``"""
    # Undocumented default values
    theory: str = "VSSTLS"


@dataclass
class Result(stls.Result):
    """
    Class used to store the results for the :obj:`qupled.vsstls.VSStls` class.
    """

    free_energy_grid: np.ndarray = None
    """Free energy grid"""
    free_energy_integrand: np.ndarray = None
    """Free energy integrand"""
    alpha: np.ndarray = None
    """Free parameter"""


@dataclass
class FreeEnergyIntegrand(hf.SerializableMixin):

    grid: np.ndarray = None
    """ Coupling parameter grid. Default = ``None``"""
    integrand: np.ndarray = None
    """ Free energy integrand. Default = ``None``"""

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


if __name__ == "__main__":
    from .mpi_worker import run_mpi_worker

    run_mpi_worker(Input, Result, VSStls.native_inputs_cls, VSStls.native_scheme_cls)
