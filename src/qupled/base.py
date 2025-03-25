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

    def to_native(self, native_input: any) -> any:
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
        for attr in self.__dict__.keys():
            if hasattr(native_scheme, attr):
                value = getattr(native_scheme, attr)
                setattr(self, attr, value) if value is not None else None

    def compute_rdf(self, rdf_grid: np.ndarray | None = None):
        if rdf_grid is None:
            self.rdf_grid = np.arange(0.0, 10.0, 0.01)
        self.rdf = native.compute_rdf(self.rdf_grid, self.wvg, self.ssf)


# -----------------------------------------------------------------------
# ClassicScheme class
# -----------------------------------------------------------------------


class ClassicScheme:

    def __init__(self):
        # File to store output on disk
        self.db_handler = DataBaseHandler()
        self.inputs: Input | None = None
        self.results: Result | None = None

    # Compute the scheme
    @util.MPI.record_time
    @util.MPI.synchronize_ranks
    def compute(self, inputs: Input, native_cls, result_cls: Result) -> None:
        self.inputs = inputs
        scheme = native_cls(self.inputs.to_native())
        status = scheme.compute()
        self._check_status_and_clean(status, scheme.recovery)
        self.results = result_cls(scheme)
        self._save()

    # Check that the dielectric scheme was solved without errors
    @util.MPI.run_only_on_root
    def _check_status_and_clean(self, status: bool, recovery: str) -> None:
        """Checks that the scheme was solved correctly and removes temporary files generated at run-time

        Args:
            status: status obtained from the native code. If status == 0 the scheme was solved correctly.
            recovery: name of the recovery file.
        """
        if status == 0:
            if os.path.isfile(recovery):
                os.remove(recovery)
            print("Dielectric theory solved successfully!")
        else:
            sys.exit("Error while solving the dielectric theory")

    @util.MPI.run_only_on_root
    def _save(self) -> None:
        if self.results is not None:
            self.db_handler.insert_run(self.inputs, self.results)

    # Compute radial distribution function
    @util.MPI.run_only_on_root
    def compute_rdf(self, rdf_grid: np.ndarray = None) -> None:
        """Computes the radial distribution function from the data stored in the output file.

        Args:
            rdf_grid: The grid used to compute the radial distribution function.
                Default = ``None`` (see :func:`qupled.util.Hdf.computeRdf`)


        """
        if self.results is not None:
            self.results.compute_rdf(rdf_grid)
            self.db_handler.insert_results_data(
                {
                    util.HDF.ResultNames.RDF.value: self.results.rdf,
                    util.HDF.ResultNames.RDF_GRID.value: self.results.rdf_grid,
                }
            )

    # Plot results
    @util.MPI.run_only_on_root
    def plot(
        self,
        to_plot: list[str],
        matsubara: np.ndarray = None,
        rdf_grid: np.ndarray = None,
    ) -> None:
        """Plots the results stored in the output file.

        Args:
            to_plot: A list of quantities to plot. This list can include all the values written to the
                 output hdf file. The radial distribution function is computed and added to the output
                 file if necessary
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Default = None, all matsubara frequencies are plotted)
            rdf_grid: The grid used to compute the radial distribution function. Applies only when the radial
                distribution function is plotted. Default = ``None`` (see :func:`qupled.util.Hdf.computeRdf`).

        """
        if util.HDF.ResultNames.RDF.value in to_plot:
            self.compute_rdf(rdf_grid)
        util.HDF.plot(
            to_plot, self.db_handler.run_id, self.db_handler.database_name, matsubara
        )


# -----------------------------------------------------------------------
# IterativeScheme class
# -----------------------------------------------------------------------


class IterativeScheme(ClassicScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(
        run_id: str, database_name: str | None = None
    ) -> IterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from a database.

        Args:
            run_id: ID of the run used to extract the information for the initial guess.
            database_name: Name of the database file. Default is None.

        Returns:
            An instance of IterativeScheme.Guess containing the initial guess data.
        """
        names = [util.HDF.ResultNames.WVG.value, util.HDF.ResultNames.SLFC.value]
        data = util.HDF.read_results(run_id, database_name, names)
        return IterativeScheme.Guess(data[names[0]], data[names[1]])

    # Initial guess
    class Guess:

        def __init__(self, wvg: np.ndarray = None, slfc: np.ndarray = None):
            self.wvg = wvg
            """ Wave-vector grid. Default = ``None``"""
            self.slfc = slfc
            """ Static local field correction. Default = ``None``"""

        def to_native(self) -> native.StlsGuess:
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
    def get_initial_guess(file_name: str) -> QuantumIterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the initial guess.
        """
        hdf_data = util.HDF.read(
            file_name,
            [
                util.HDF.ResultNames.WVG.value,
                util.HDF.ResultNames.SSF.value,
                util.HDF.ResultNames.ADR.value,
                util.HDF.ResultNames.MATSUBARA.value,
            ],
        )
        return QuantumIterativeScheme.Guess(
            hdf_data[util.HDF.ResultNames.WVG.value],
            hdf_data[util.HDF.ResultNames.SSF.value],
            np.ascontiguousarray(hdf_data[util.HDF.ResultNames.ADR.value]),
            hdf_data[util.HDF.ResultNames.MATSUBARA.value],
        )

    # Save results to disk
    @util.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        if scheme.inputs.degeneracy > 0:
            pd.DataFrame(scheme.adr).to_hdf(
                self.hdf_file_name, key=util.HDF.ResultNames.ADR.value
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
            native_guess = native.QstlsGuess()
            for attr, value in self.__dict__.items():
                native_value = value if value is not None else np.empty(0)
                setattr(native_guess, attr, native_value)
            return native_guess
