from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

from . import native
from . import util
from qupled.database import DataBaseHandler


# -----------------------------------------------------------------------
# Input class
# -----------------------------------------------------------------------


class Input:
    """
    Class representing the input parameters for a computational model.
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

    def to_native(self, native_input: any) -> any:
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

        Returns:
            any: The `native_input` object with updated attributes.
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
        return native_input


# -----------------------------------------------------------------------
# Result class
# -----------------------------------------------------------------------


class Result:
    """
    Class representing the results for a computational model.
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
        self.rdf_grid = rdf_grid if rdf_grid is not None else np.arange(0.0, 10.0, 0.01)
        self.rdf = native.compute_rdf(self.rdf_grid, self.wvg, self.ssf)


# -----------------------------------------------------------------------
# ClassicScheme class
# -----------------------------------------------------------------------


class ClassicScheme:

    def __init__(self):
        # File to store output on disk
        self.db_handler = DataBaseHandler()
        self.inputs: Input = None
        """An object representing the inputs. Default = ``None``"""
        self.results: Result = None
        """An object representing the results. Default = ``None``"""

    @property
    def run_id(self):
        """
        Property that retrieves the run ID from the database handler.

        Returns:
            str: The run ID associated with the current database handler.
        """
        return self.db_handler.run_id

    # Compute the scheme
    @util.MPI.record_time
    @util.MPI.synchronize_ranks
    def compute(
        self, inputs: Input, native_scheme_cls: any, native_input: any, result: Result
    ):
        """
        Perform computation using the provided inputs, scheme, and result objects.

        Args:
            inputs (Input): The input data for the computation.
            native_scheme_cls (any): The class representing the native computation scheme.
            native_input (any): The native input data required by the computation scheme.
            result (Result): The object to store the computation results.

        Notes:
            - Converts the provided inputs to a native format compatible with the computation scheme.
            - Initializes the computation scheme with the native inputs and executes the computation.
            - Checks the computation status and performs recovery if necessary.
            - Populates the result object with the computation results and saves the state.
        """
        self.inputs = inputs
        native_inputs = self.inputs.to_native(native_input)
        scheme = native_scheme_cls(native_inputs)
        status = scheme.compute()
        self._check_status_and_clean(status, scheme.recovery)
        self.results = result
        self.results.from_native(scheme)
        self._save()

    # Check that the dielectric scheme was solved without errors
    @util.MPI.run_only_on_root
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

    @util.MPI.run_only_on_root
    def _save(self):
        """
        Saves the current run's inputs and results to the database.

        This method checks if the `results` attribute is not None, and if so,
        it uses the `db_handler` to insert the current run's `inputs` and `results`
        into the database.
        """
        if self.results is not None:
            self.db_handler.insert_run(self.inputs, self.results)

    # Compute radial distribution function
    @util.MPI.run_only_on_root
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


# -----------------------------------------------------------------------
# IterativeScheme class
# -----------------------------------------------------------------------


class IterativeScheme(ClassicScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(run_id: str, database_name: str | None = None) -> Guess:
        """Constructs an initial guess object by extracting the information from a database.

        Args:
            run_id: ID of the run used to extract the information for the initial guess.
            database_name: Name of the database file. Default is None.

        Returns:
            An instance of Guess containing the initial guess data.
        """
        names = ["wvg", "slfc"]
        data = util.DataBase.read_results(run_id, database_name, names)
        return IterativeScheme.Guess(data[names[0]], data[names[1]])

    # Initial guess
    class Guess:

        def __init__(self, wvg: np.ndarray = None, slfc: np.ndarray = None):
            self.wvg = wvg
            """ Wave-vector grid. Default = ``None``"""
            self.slfc = slfc
            """ Static local field correction. Default = ``None``"""

        def to_native(self) -> native.StlsGuess:
            """
            Converts the current object to a native `StlsGuess` object.

            This method iterates over the attributes of the current object and
            assigns their values to a new `StlsGuess` object. If an attribute's
            value is `None`, it is replaced with an empty NumPy array.

            Returns:
                native.StlsGuess: A new instance of `StlsGuess` with attributes
                copied from the current object.
            """
            native_guess = native.StlsGuess()
            for attr, value in self.__dict__.items():
                native_value = value if value is not None else np.empty(0)
                setattr(native_guess, attr, native_value)
            return native_guess


# -----------------------------------------------------------------------
# QuantumIterativeScheme class
# -----------------------------------------------------------------------


class QuantumIterativeScheme(IterativeScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(run_id: str, database_name: str | None = None) -> Guess:
        """
        Generates an initial guess for a quantum iterative scheme based on data
        retrieved from a database.

        Args:
            run_id (str): The unique identifier for the run whose data is to be retrieved.
            database_name (str | None, optional): The name of the database to query.
                If None, a default database is used.

        Returns:
            QuantumIterativeScheme.Guess: An object containing the initial guess values
            for the quantum iterative scheme, including results and inputs.

        Raises:
            Exception: If there is an issue reading data from the database.
        """
        result_names = ["wvg", "ssf", "adr"]
        input_names = ["matsubara"]
        data = util.DataBase.read_run(run_id, database_name, input_names, result_names)
        inputs = data[DataBaseHandler.INPUTS_TABLE_NAME]
        results = data[DataBaseHandler.RESULTS_TABLE_NAME]
        return QuantumIterativeScheme.Guess(
            results[result_names[0]],
            results[result_names[1]],
            results[result_names[2]],
            inputs[input_names[0]],
        )

    # Initial guess
    class Guess:

        def __init__(
            self,
            wvg: np.ndarray = None,
            ssf: np.ndarray = None,
            adr: np.ndarray = None,
            matsubara: int = 0,
        ):
            self.wvg = wvg
            """ Wave-vector grid. Default = ``None``"""
            self.ssf = ssf
            """ Static structure factor. Default = ``None``"""
            self.adr = adr
            """ Auxiliary density response. Default = ``None``"""
            self.matsubara = matsubara
            """ Number of matsubara frequencies. Default = ``0``"""

        def to_native(self) -> native.QStlsGuess:
            """
            Converts the current object to a native `QStlsGuess` object.

            This method creates an instance of `native.QStlsGuess` and populates its
            attributes with the corresponding values from the current object's
            attributes. If an attribute's value is `None`, it is replaced with an
            empty NumPy array.

            Returns:
                native.QStlsGuess: A new instance of `native.QStlsGuess` with attributes
                populated from the current object.
            """
            native_guess = native.QstlsGuess()
            for attr, value in self.__dict__.items():
                native_value = value if value is not None else np.empty(0)
                setattr(native_guess, attr, native_value)
            return native_guess
