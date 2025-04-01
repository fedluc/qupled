from __future__ import annotations
import os
import numpy as np
from . import native
from . import mpi
from . import database


class Rpa:
    """
    Class used to solve the RPA scheme.
    """

    def __init__(self):
        self.inputs: Input = None
        """The inputs used to solve the scheme. Default = ``None``"""
        self.results: Result = Result()
        """The results obtained by solving the scheme"""
        # Undocumented properties
        self.db_handler = database.DataBaseHandler()
        self.native_scheme_cls = native.Rpa
        self.native_inputs = native.RpaInput()

    @property
    def run_id(self):
        """
        Property that retrieves the run ID from the database handler.

        Returns:
            str: The run ID associated with the current database handler.
        """
        return self.db_handler.run_id

    # Compute
    def compute(self, inputs: Input):
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self.inputs = inputs
        self.inputs.to_native(self.native_inputs)
        scheme = self.native_scheme_cls(self.native_inputs)
        status = scheme.compute()
        self._check_status_and_clean(status, scheme.recovery)
        self.results.from_native(scheme)
        self._save()

    # Compute radial distribution function
    @mpi.MPI.run_only_on_root
    def compute_rdf(self, rdf_grid: np.ndarray = None):
        """
        Computes the radial distribution function (RDF) using the provided RDF grid.
        If results are available, this method computes the RDF and stores the results
        in the database.

        Args:
            rdf_grid: A numpy array representing the RDF grid.
                If not provided, a default grid will be used.
        """
        if self.results is not None:
            self.results.compute_rdf(rdf_grid)
            self.db_handler.insert_results(
                {"rdf": self.results.rdf, "rdf_grid": self.results.rdf_grid}
            )

    # Check that the dielectric scheme was solved without errors
    @mpi.MPI.run_only_on_root
    def _check_status_and_clean(self, status: bool, recovery: str):
        """
        Checks the status of a process and performs cleanup if successful.

        Args:
            status (bool): The status of the process. A value of 0 indicates success.
            recovery (str): The file path to a recovery file that should be removed
                            if the process is successful.

        Raises:
            RuntimeError: If the status is not 0, indicating an error in the process.

        Side Effects:
            - Removes the recovery file if it exists and the status is 0.
            - Prints a success message if the status is 0.
        """
        if status == 0:
            if os.path.exists(recovery):
                os.remove(recovery)
            print("Dielectric theory solved successfully!")
        else:
            raise RuntimeError("Error while solving the dielectric theory")

    @mpi.MPI.run_only_on_root
    def _save(self):
        """
        Saves the current run's inputs and results to the database.

        This method checks if the `results` attribute is not None, and if so,
        it uses the `db_handler` to insert the current run's `inputs` and `results`
        into the database.
        """
        self.db_handler.insert_run(self.inputs, self.results)


class Input:
    """
    Class used to store the inputs for the :obj:`qupled.rpa.Rpa` class.
    """

    def __init__(self, coupling: float, degeneracy: float):
        """
        Initialize the base class with the given parameters.

        Parameters:
            coupling (float): Coupling parameter.
            degeneracy (float): Degeneracy parameter.
        """
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

        - full: the inner integral is evaluated at arbitrary points selected automatically by the quadrature rule

        - segregated: the inner integral is evaluated on a fixed grid that depends on the integrand that is being processed

        Segregated is usually faster than full but it could become
        less accurate if the fixed points are not chosen correctly. Default =  ``'full'``
        """
        self.matsubara: int = 128
        """Number of Matsubara frequencies. Default = ``128``"""
        self.resolution: float = 0.1
        """Resolution of the wave-vector grid. Default =  ``0.1``"""
        self.threads: int = 1
        """Number of OMP threads for parallel calculations. Default =  ``1``"""
        # Undocumented default values
        self.theory: str = "RPA"

    def to_native(self, native_input: any):
        """
        Converts the attributes of the current object to their native representations
        and sets them on the provided `native_input` object.

        This method iterates through the attributes of the current object and checks
        if the `native_input` object has a corresponding attribute. If it does, the
        method attempts to convert the attribute's value to its native representation
        using a `to_native` method, if available. Otherwise, it directly assigns the
        attribute's value to the `native_input` object.

        Args:
            native_input (any): The object to which the native representations of the
                current object's attributes will be assigned.
        """
        name = Input.to_native.__name__
        for attr, value in self.__dict__.items():
            if hasattr(native_input, attr):
                value_to_set = (
                    tonative()
                    if callable(tonative := getattr(value, name, None))
                    else value
                )
                setattr(native_input, attr, value_to_set)


class Result:
    """
    Class used to store the results for the :obj:`qupled.rpa.Rpa` class.
    """

    def __init__(self):
        self.idr: np.ndarray = None
        """Ideal density response"""
        self.rdf: np.ndarray = None
        """Radial distribution function"""
        self.rdf_grid: np.ndarray = None
        """Radial distribution function grid"""
        self.sdr: np.ndarray = None
        """Static density response"""
        self.slfc: np.ndarray = None
        """Static local field correction"""
        self.ssf: np.ndarray = None
        """Static structure factor"""
        self.uint: float = None
        """Internal energy"""
        self.wvg: np.ndarray = None
        """Wave-vector grid"""

    def from_native(self, native_scheme: any):
        """
        Updates the attributes of the current object based on the attributes of a given native scheme object.

        Args:
            native_scheme (any): An object containing attributes to update the current object with.

        Notes:
            - Only attributes that exist in both the current object and the native_scheme object will be updated.
            - Attributes with a value of `None` in the native_scheme object will not overwrite the current object's attributes.
        """
        for attr in self.__dict__.keys():
            if hasattr(native_scheme, attr):
                value = getattr(native_scheme, attr)
                setattr(self, attr, value) if value is not None else None

    def compute_rdf(self, rdf_grid: np.ndarray | None = None):
        """
        Compute the radial distribution function (RDF) for the system.

        Args:
            rdf_grid (np.ndarray | None, optional): A 1D array specifying the grid points
                at which the RDF is computed. If None, a default grid ranging from 0.0
                to 10.0 with a step size of 0.01 is used.

        Returns:
            None: The computed RDF is stored in the `self.rdf` attribute.
        """
        if self.wvg is not None and self.ssf is not None:
            self.rdf_grid = (
                rdf_grid if rdf_grid is not None else np.arange(0.0, 10.0, 0.01)
            )
            self.rdf = native.compute_rdf(self.rdf_grid, self.wvg, self.ssf)
