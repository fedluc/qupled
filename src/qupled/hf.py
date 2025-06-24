from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np

from . import database
from . import mpi
from . import native


class HF:
    """
    Class used to solve the HF scheme.
    """

    # Mapping of native scheme status to run status in the database
    NATIVE_TO_RUN_STATUS = {
        0: database.DataBaseHandler.RunStatus.SUCCESS,
        1: database.DataBaseHandler.RunStatus.FAILED,
    }

    # Native classes used to solve the scheme
    native_scheme_cls = native.HF
    native_inputs_cls = native.Input

    def __init__(self):
        self.inputs: Input = None
        """The inputs used to solve the scheme. Default = ``None``"""
        self.results: Result = Result()
        """The results obtained by solving the scheme"""
        # Undocumented properties
        self.db_handler = database.DataBaseHandler()
        self.native_scheme_status = None

    @property
    def run_id(self):
        """
        Property that retrieves the run ID from the database handler.

        Returns:
            str: The run ID associated with the current database handler.
        """
        return self.db_handler.run_id

    @mpi.MPI.record_time
    @mpi.MPI.synchronize_ranks
    def compute(self, inputs: Input):
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self.inputs = inputs
        self._add_run_to_database()
        self._compute_native_new()
        self._save()

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
                {"rdf": self.results.rdf, "rdf_grid": self.results.rdf_grid},
                conflict_mode=database.DataBaseHandler.ConflictMode.UPDATE,
            )

    def _add_run_to_database(self):
        """
        Adds the current run information to the database.

        This method inserts the run details stored in `self.inputs` into the database
        using the `db_handler`. It also updates the `database_info` attribute of
        `self.inputs` with the current `run_id`.
        """
        self.db_handler.insert_run(self.inputs)
        self.inputs.database_info.run_id = self.run_id

    def _compute_native(self):
        """
        Computes the native representation of the inputs and processes the results.

        This method performs the following steps:
        1. Converts the current inputs to their native representation.
        2. Initializes a native scheme object using the native inputs.
        3. Computes the native scheme and stores its status.
        4. Converts the results from the native scheme back to the desired format.
        """
        native_inputs = self.native_inputs_cls()
        self.inputs.to_native(native_inputs)
        scheme = self.native_scheme_cls(native_inputs)
        self.native_scheme_status = scheme.compute()
        self.results.from_native(scheme)

    def _compute_native_new(self):
        input_dict = self.inputs.to_dict()
        input_path = Path("input.json")
        result_path = Path("results.json")
        with input_path.open("w") as f:
            json.dump(input_dict, f)
        subprocess.run(
            ["mpiexec", "-n", "1", "python", "-m", self.__module__], check=True
        )
        with result_path.open() as f:
            result_dict = json.load(f)
        self.results = Result.from_dict(result_dict)

    @mpi.MPI.run_only_on_root
    def _save(self):
        """
        Saves the current state and results to the database.

        This method updates the run status in the database using the current
        native scheme status and inserts the results into the database.
        """
        run_status = self.NATIVE_TO_RUN_STATUS.get(
            self.native_scheme_status, database.DataBaseHandler.RunStatus.FAILED
        )
        self.db_handler.update_run_status(run_status)
        self.db_handler.insert_results(self.results.__dict__)


class SerializableMixin:
    def to_dict(self):
        result = {}
        for key, value in self.__dict__.items():
            if hasattr(value, "to_dict") and callable(value.to_dict):
                result[key] = value.to_dict()
            elif isinstance(value, np.ndarray):
                result[key] = value.tolist()
            else:
                result[key] = value
        return result


class Input(SerializableMixin):
    """
    Class used to store the inputs for the :obj:`qupled.hf.HF` class.
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
        self.theory: str = "HF"
        self.database_info: DatabaseInfo = DatabaseInfo()

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
            if hasattr(native_input, attr) and value is not None:
                value_to_set = (
                    tonative()
                    if callable(tonative := getattr(value, name, None))
                    else value
                )
                setattr(native_input, attr, value_to_set)

    @classmethod
    def from_dict(cls, d):
        obj = cls.__new__(cls)
        for key, value in d.items():
            if key == "database_info" and isinstance(value, dict):
                obj.database_info = DatabaseInfo.from_dict(value)
            else:
                setattr(obj, key, value)

        return obj


class Result(SerializableMixin):
    """
    Class used to store the results for the :obj:`qupled.hf.HF` class.
    """

    def __init__(self):
        self.idr: np.ndarray = None
        """Ideal density response"""
        self.lfc: np.ndarray = None
        """Local field correction"""
        self.rdf: np.ndarray = None
        """Radial distribution function"""
        # self.rdf_grid: np.ndarray = None
        # """Radial distribution function grid"""
        self.sdr: np.ndarray = None
        """Static density response"""
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
                valid_value = value is not None and not callable(value)
                setattr(self, attr, value) if valid_value else None

    @classmethod
    def from_dict(cls, d):
        obj = cls.__new__(cls)
        for key, value in d.items():
            if isinstance(value, list):  # assume these were NumPy arrays
                setattr(obj, key, np.array(value))
            else:
                setattr(obj, key, value)
        return obj

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


class DatabaseInfo(SerializableMixin):
    """
    Class used to store the database information passed to the native code.
    """

    def __init__(self):
        self.name: str = database.DataBaseHandler.DEFAULT_DATABASE_NAME
        """Database name"""
        self.run_id: int = None
        """ID of the run in the database"""
        self.run_table_name: str = database.DataBaseHandler.RUN_TABLE_NAME
        """Name of the table used to store the runs in the database"""

    def to_native(self) -> native.DatabaseInfo:
        """
        Converts the current object to a native `DatabaseInfo` instance.
        This method creates a new instance of `native.DatabaseInfo` and copies
        all non-None attributes from the current object to the new instance.
        Returns:
            native.DatabaseInfo: A new instance of `native.DatabaseInfo` with
            attributes copied from the current object.
        """
        native_database_info = native.DatabaseInfo()
        for attr, value in self.__dict__.items():
            if value is not None:
                setattr(native_database_info, attr, value)
        return native_database_info

    @classmethod
    def from_dict(cls, d):
        obj = cls.__new__(cls)
        for key, value in d.items():
            setattr(obj, key, value)
        return obj


if __name__ == "__main__":
    from .mpi_worker import run_mpi_worker

    run_mpi_worker(Input, Result, HF.native_inputs_cls, HF.native_scheme_cls)
