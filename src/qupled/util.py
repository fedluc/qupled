import functools
from enum import Enum

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps as cm

from . import native
from . import database

# -----------------------------------------------------------------------
# Hdf class
# -----------------------------------------------------------------------


class HDF:
    """Class to manipulate the output hdf files produced when a scheme is solved."""

    class EntryKeys(Enum):
        ALPHA = "alpha"
        ADR = "adr"
        BF = "bf"
        COUPLING = "coupling"
        CUTOFF = "cutoff"
        FREQUENCY_CUTOFF = "frequency_cutoff"
        DEGENERACY = "degeneracy"
        ERROR = "error"
        FXC_GRID = "fxc_grid"
        FXC_INT = "fxc_int"
        MATSUBARA = "matsubara"
        IDR = "idr"
        INFO = "info"
        RESOLUTION = "resolution"
        RDF = "rdf"
        RDF_GRID = "rdf_grid"
        SDR = "sdr"
        SLFC = "slfc"
        SSF = "ssf"
        SSF_HF = "ssf_HF"
        THEORY = "theory"
        WVG = "wvg"

    class EntryType(Enum):
        NUMPY = "numpy"
        NUMPY2D = "numpy2D"
        NUMBER = "number"
        STRING = "string"

    # Construct
    def __init__(self):
        self.entries = {
            HDF.EntryKeys.ALPHA.value: self.Entries(
                "Free Parameter for VS schemes", HDF.EntryType.NUMPY.value
            ),
            HDF.EntryKeys.ADR.value: self.Entries(
                "Auxiliary density response", HDF.EntryType.NUMPY2D.value
            ),
            HDF.EntryKeys.BF.value: self.Entries(
                "Bridge function adder", HDF.EntryType.NUMPY.value
            ),
            HDF.EntryKeys.COUPLING.value: self.Entries(
                "Coupling parameter", HDF.EntryType.NUMBER.value
            ),
            HDF.EntryKeys.CUTOFF.value: self.Entries(
                "Cutoff for the wave-vector grid", HDF.EntryType.NUMBER.value
            ),
            HDF.EntryKeys.FREQUENCY_CUTOFF.value: self.Entries(
                "Cutoff for the frequency", HDF.EntryType.NUMBER.value
            ),
            HDF.EntryKeys.DEGENERACY.value: self.Entries(
                "Degeneracy parameter", HDF.EntryType.NUMBER.value
            ),
            HDF.EntryKeys.ERROR.value: self.Entries(
                "Residual error in the solution", HDF.EntryType.NUMBER.value
            ),
            HDF.EntryKeys.FXC_GRID.value: self.Entries(
                "Coupling parameter", HDF.EntryType.NUMPY.value
            ),
            HDF.EntryKeys.FXC_INT.value: self.Entries(
                "Free Energy integrand", HDF.EntryType.NUMPY2D.value
            ),
            HDF.EntryKeys.MATSUBARA.value: self.Entries(
                "Number of matsubara frequencies", HDF.EntryType.NUMBER.value
            ),
            HDF.EntryKeys.IDR.value: self.Entries(
                "Ideal density response", HDF.EntryType.NUMPY2D.value
            ),
            HDF.EntryKeys.RESOLUTION.value: self.Entries(
                "Resolution for the wave-vector grid", HDF.EntryType.NUMBER.value
            ),
            HDF.EntryKeys.RDF.value: self.Entries(
                "Radial distribution function", HDF.EntryType.NUMPY.value
            ),
            HDF.EntryKeys.RDF_GRID.value: self.Entries(
                "Inter-particle distance", HDF.EntryType.NUMPY.value
            ),
            HDF.EntryKeys.SDR.value: self.Entries(
                "Static density response", HDF.EntryType.NUMPY.value
            ),
            HDF.EntryKeys.SLFC.value: self.Entries(
                "Static local field correction", HDF.EntryType.NUMPY.value
            ),
            HDF.EntryKeys.SSF.value: self.Entries(
                "Static structure factor", HDF.EntryType.NUMPY.value
            ),
            HDF.EntryKeys.SSF_HF.value: self.Entries(
                "Hartree-Fock static structure factor", HDF.EntryType.NUMPY.value
            ),
            HDF.EntryKeys.THEORY.value: self.Entries(
                "Theory that is being solved", HDF.EntryType.STRING.value
            ),
            HDF.EntryKeys.WVG.value: self.Entries(
                "Wave-vector", HDF.EntryType.NUMPY.value
            ),
        }

    # Structure used to categorize the entries stored in the hdf file
    class Entries:
        def __init__(self, description, entry_type):
            self.description = description  # Descriptive string of the entry
            self.entry_type = (
                entry_type  # Type of entry (numpy, numpy2, number or string)
            )
            assert (
                self.entry_type == HDF.EntryType.NUMPY.value
                or self.entry_type == HDF.EntryType.NUMPY2D.value
                or self.entry_type == HDF.EntryType.NUMBER.value
                or self.entry_type == HDF.EntryType.STRING.value
            )

    # Read runs in the database
    def read_runs(self, database_name: str | None = None) -> dict:
        """Reads runs from the database and returns the content in the form of a dictionary.

        Args:
            database_name (str, optional): Name of the database to read from. Defaults to None.

        Returns:
            A dictionary whose keys are the run ids and values are the corresponding runs information.
        """
        db_handler = database.DataBaseHandler(database_name)
        return db_handler.get_runs()

    # Read inputs in the database
    def read_inputs(
        self,
        run_id: int,
        database_name: str | None = None,
        names: list[str] | None = None,
    ) -> dict:
        """Reads inputs from the database and returns the content in the form of a dictionary.

        Args:
            run_id: Identifier of the run to read input for.
            database_name: Name of the database to read from (default is None).
            names: A list of quantities to read (default is None, which reads all available quantities).

        Returns:
            A dictionary whose keys are the quantities listed in names and values are the corresponding inputs.
        """
        db_handler = database.DataBaseHandler(database_name)
        return db_handler.get_inputs_data(run_id, names if names is not None else [])

    # Read results in the database
    def read_results(
        self,
        run_id: int,
        database_name: str | None = None,
        names: list[str] | None = None,
    ) -> dict:
        """Reads results from the database and returns the content in the form of a dictionary.

        Args:
            run_id: Identifier of the run to read results for.
            database_name: Name of the database to read from (default is None).
            names: A list of quantities to read (default is None, which reads all available quantities).

        Returns:
            A dictionary whose keys are the quantities listed in names and values are the corresponding results.
        """
        db_handler = database.DataBaseHandler(database_name)
        return db_handler.get_results_data(run_id, names)

    # Plot from data in hdf file
    def plot(
        self,
        to_plot: list[str],
        run_id: id,
        database_name: str | None = None,
        matsubara: np.array = None,
    ) -> None:
        """Plots the results stored in an hdf file.

        Args:
            hdf: Name of the hdf file
            to_plot: A list of quantities to plot. Allowed quantities include adr (auxiliary density response),
                bf (bridge function adder), fxci (free energy integrand), idr (ideal density response), rdf
                (radial distribution function), sdr (static density response), slfc (static local field correction)
                ssf (static structure factor) and ssfHF (Hartree-Fock static structure factor).
                If the hdf file does not contain the specified quantity, an error is thrown
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Defaults to  None, all matsubara frequencies are plotted)

        """
        for name in to_plot:
            description = (
                self.entries[name].description if name in self.entries.keys() else ""
            )
            if name == HDF.EntryKeys.RDF.value:
                x = self.read_results(
                    run_id, database_name, [name, HDF.EntryKeys.RDF_GRID.value]
                )
                Plot.plot_1D(
                    x[HDF.EntryKeys.RDF_GRID.value],
                    x[name],
                    self.entries[HDF.EntryKeys.RDF_GRID.value].description,
                    description,
                )
            elif name in [HDF.EntryKeys.ADR.value, HDF.EntryKeys.IDR.value]:
                x = self.read_results(
                    run_id,
                    database_name,
                    [name, HDF.EntryKeys.WVG.value, HDF.EntryKeys.MATSUBARA.value],
                )
                if matsubara is None:
                    matsubara = np.arange(x[HDF.EntryKeys.MATSUBARA.value])
                Plot.plot_1D_parametric(
                    x[HDF.EntryKeys.WVG.value],
                    x[name],
                    self.entries[HDF.EntryKeys.WVG.value].description,
                    description,
                    matsubara,
                )
            elif name == HDF.EntryKeys.FXC_INT.value:
                x = self.read_results(
                    run_id, database_name, [name, HDF.EntryKeys.FXC_GRID.value]
                )
                Plot.plot_1D(
                    x[HDF.EntryKeys.FXC_GRID.value],
                    x[name][1, :],
                    self.entries[HDF.EntryKeys.FXC_GRID.value].description,
                    description,
                )
            elif name in [
                HDF.EntryKeys.BF.value,
                HDF.EntryKeys.SDR.value,
                HDF.EntryKeys.SLFC.value,
                HDF.EntryKeys.SSF.value,
                HDF.EntryKeys.SSF_HF.value,
            ]:
                x = self.read_results(
                    run_id, database_name, [name, HDF.EntryKeys.WVG.value]
                )
                Plot.plot_1D(
                    x[HDF.EntryKeys.WVG.value],
                    x[name],
                    self.entries[HDF.EntryKeys.WVG.value].description,
                    self.entries[name].description,
                )
            elif name == HDF.EntryKeys.ALPHA.value:
                x = self.read_results(
                    run_id, database_name, [name, HDF.EntryKeys.FXC_GRID.value]
                )
                Plot.plot_1D(
                    x[HDF.EntryKeys.FXC_GRID.value][::2],
                    x[name][::2],
                    self.entries[HDF.EntryKeys.FXC_GRID.value].description,
                    self.entries[name].description,
                )
            else:
                raise ValueError(f"Unknown quantity to plot: {name}")


# -----------------------------------------------------------------------
# Plot class
# -----------------------------------------------------------------------


class Plot:
    """Class to collect methods used for plotting"""

    # One dimensional plots
    @staticmethod
    def plot_1D(x, y, xlabel, ylabel):
        """Produces the plot of a one-dimensional quantity.

        Positional arguments:
        x -- data for the x-axis (a numpy array)
        y -- data for the y-axis (a numpy array)
        xlabel -- label for the x-axis (a string)
        ylabel -- label for the y-axis (a string)
        """
        plt.plot(x, y, "b")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()

    # One dimensional plots with one parameter
    @staticmethod
    def plot_1D_parametric(x, y, xlabel, ylabel, parameters):
        """Produces the plot of a one-dimensional quantity that depends on an external parameter.

        Positional arguments:
        x -- data for the x-axis (a numpy array)
        y -- data for the y-axis (a two-dimensional numpy array)
        xlabel -- label for the x-axis (a string)
        ylabel -- label for the y-axis (a string)
        parameters -- list of parameters for which the results should be plotted
        """
        num_parameters = parameters.size
        cmap = cm["viridis"]
        for i in np.arange(num_parameters):
            color = cmap(1.0 * i / num_parameters)
            plt.plot(x, y[:, parameters[i]], color=color)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()


# -----------------------------------------------------------------------
# MPI class
# -----------------------------------------------------------------------


class MPI:
    """Class to handle the calls to the MPI API"""

    def __init__(self):
        self.qp_mpi = native.MPI

    def rank(self):
        """Get rank of the process"""
        return self.qp_mpi.rank()

    def is_root(self):
        """Check if the current process is root (rank 0)"""
        return self.qp_mpi.is_root()

    def barrier(self):
        """Setup an MPI barrier"""
        self.qp_mpi.barrier()

    def timer(self):
        """Get wall time"""
        return self.qp_mpi.timer()

    @staticmethod
    def run_only_on_root(func):
        """Python decorator for all methods that have to be run only by root"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if MPI().is_root():
                return func(*args, **kwargs)

        return wrapper

    @staticmethod
    def synchronize_ranks(func):
        """Python decorator for all methods that need rank synchronization"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            func(*args, **kwargs)
            MPI().barrier()

        return wrapper

    @staticmethod
    def record_time(func):
        """Python decorator for all methods that have to be timed"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            mpi = MPI()
            tic = mpi.timer()
            func(*args, **kwargs)
            toc = mpi.timer()
            dt = toc - tic
            hours = dt // 3600
            minutes = (dt % 3600) // 60
            seconds = dt % 60
            if mpi.is_root():
                if hours > 0:
                    print("Elapsed time: %d h, %d m, %d s." % (hours, minutes, seconds))
                elif minutes > 0:
                    print("Elapsed time: %d m, %d s." % (minutes, seconds))
                else:
                    print("Elapsed time: %.1f s." % seconds)

        return wrapper
