from __future__ import annotations

from dataclasses import field

import numpy as np

from . import database
from . import mpi
from . import native
from . import serialize
from . import timer


class Solver:
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

    @timer.timer
    def compute(self, inputs: Input):
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self.inputs = inputs
        self._add_run_to_database()
        self._compute_native()
        self._save()

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
            self.results.compute_rdf(self.inputs.dimension, rdf_grid)
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
        self.inputs.database_info = DatabaseInfo(
            blob_storage=self.db_handler.blob_storage,
            name=self.db_handler.engine.url.database,
            run_id=self.run_id,
        )

    def _compute_native(self):
        """
        Determines whether to execute the native computation in parallel or serial mode.

        Checks if MPI (Message Passing Interface) is available and if the number of requested processes
        is greater than one. If both conditions are met, runs the computation in parallel; otherwise,
        runs it in serial mode.
        """
        if native.uses_mpi:
            self._compute_native_mpi()
        else:
            self._compute_native_serial()

    def _compute_native_serial(self):
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

    def _compute_native_mpi(self):
        """
        Executes a native MPI computation workflow.

        This method performs the following steps:
        1. Writes the necessary input files for the MPI computation using `mpi.write_inputs`.
        2. Launches the MPI execution by calling `mpi.launch_mpi_execution` with the current module and the specified number of processes.
        3. Reads the computation results using `mpi.read_results` and assigns them to `self.results`.
        4. Cleans up any temporary files generated during the computation with `mpi.clean_files`.
        """
        mpi.write_inputs(self.inputs)
        mpi.launch_mpi_execution(self.__module__, self.inputs.processes)
        self.native_scheme_status = mpi.read_status()
        self.results = mpi.read_results(type(self.results))
        mpi.clean_files()

    @classmethod
    def run_mpi_worker(cls, InputCls, ResultCls):
        inputs = mpi.read_inputs(InputCls)
        native_inputs = cls.native_inputs_cls()
        inputs.to_native(native_inputs)
        scheme = cls.native_scheme_cls(native_inputs)
        status = scheme.compute()
        mpi.write_results(scheme, ResultCls)
        mpi.write_status(scheme, status)

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


@serialize.serializable_dataclass
class Input:
    """
    Class used to store the inputs for the :obj:`qupled.hf.HF` class.
    """

    coupling: float
    """Coupling parameter."""
    degeneracy: float
    """Degeneracy parameter."""
    chemical_potential: list[float] = field(default_factory=lambda: [-10.0, 10.0])
    """Initial guess for the chemical potential. Default = ``[-10, 10]``"""
    cutoff: float = 10.0
    """Cutoff for the wave-vector grid. Default =  ``10.0``"""
    dimension: str = "D3"
    """Dimesion of the system."""
    frequency_cutoff: float = 10.0
    """Cutoff for the frequency (applies only in the ground state). Default =  ``10.0``"""
    integral_error: float = 1.0e-5
    """Accuracy (relative error) in the computation of integrals. Default = ``1.0e-5``"""
    integral_strategy: str = "full"
    """
    Scheme used to solve two-dimensional integrals
    allowed options include:

    - full: the inner integral is evaluated at arbitrary points selected automatically by the quadrature rule

    - segregated: the inner integral is evaluated on a fixed grid that depends on the integrand that is being processed

    Segregated is usually faster than full but it could become
    less accurate if the fixed points are not chosen correctly. Default =  ``'full'``
    """
    matsubara: int = 128
    """Number of Matsubara frequencies. Default = ``128``"""
    resolution: float = 0.1
    """Resolution of the wave-vector grid. Default =  ``0.1``"""
    threads: int = 1
    """Number of OMP threads for parallel calculations. Default =  ``1``"""
    processes: int = 1
    """Number of MPI processes for parallel calculations. Default =  ``1``"""
    theory: str = "HF"
    database_info: DatabaseInfo = None

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
                if callable(tonative := getattr(value, name, None)):
                    value_to_set = tonative()
                elif attr == "dimension":
                    value_to_set = getattr(native.Dimension, value)
                else:
                    value_to_set = value
                setattr(native_input, attr, value_to_set)


@serialize.serializable_dataclass
class Result:
    """
    Class used to store the results for the :obj:`qupled.hf.HF` class.
    """

    idr: np.ndarray = None
    """Ideal density response"""
    lfc: np.ndarray = None
    """Local field correction"""
    rdf: np.ndarray = None
    """Radial distribution function"""
    rdf_grid: np.ndarray = None
    """Radial distribution function grid"""
    sdr: np.ndarray = None
    """Static density response"""
    ssf: np.ndarray = None
    """Static structure factor"""
    uint: float = None
    """Internal energy"""
    wvg: np.ndarray = None
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
        for attr in self.__dataclass_fields__:
            if hasattr(native_scheme, attr):
                value = getattr(native_scheme, attr)
                valid_value = value is not None and not callable(value)
                setattr(self, attr, value) if valid_value else None

    def compute_rdf(self, dimension: str, rdf_grid: np.ndarray | None = None):
        """
        Compute the radial distribution function (RDF) for the system.

        Args:
            rdf_grid (np.ndarray | None, optional): A 1D array specifying the grid points
                at which the RDF is computed. If None, a default grid ranging from 0.0
                to 10.0 with a step size of 0.01 is used.

        Returns:
            None: The computed RDF is stored in the `self.rdf` attribute.
        """
        native_dimension = getattr(native.Dimension, dimension)
        if self.wvg is not None and self.ssf is not None:
            self.rdf_grid = (
                rdf_grid if rdf_grid is not None else np.arange(0.0, 10.0, 0.01)
            )
            self.rdf = native.compute_rdf(
                self.rdf_grid, self.wvg, self.ssf, native_dimension
            )


@serialize.serializable_dataclass
class DatabaseInfo:
    """
    Class used to store the database information passed to the native code.
    """

    blob_storage: str = None
    """Directory used to store the blob data"""
    name: str = None
    """Database name"""
    run_id: int = None
    """ID of the run in the database"""
    run_table_name: str = database.DataBaseHandler.RUN_TABLE_NAME
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


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, Result)
