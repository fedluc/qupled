from __future__ import annotations

import numpy as np

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
        grid: np.ndarray = None,
        integrand: np.ndarray = None,
        alpha: np.ndarray = None,
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
