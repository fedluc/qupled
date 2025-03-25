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

    class ResultNames(Enum):
        ALPHA = "alpha"
        ADR = "adr"
        BF = "bf"
        FXC_GRID = "free_energy_grid"
        FXC_INT = "free_energy_integrand"
        IDR = "idr"
        RDF = "rdf"
        RDF_GRID = "rdf_grid"
        SDR = "sdr"
        SLFC = "slfc"
        SSF = "ssf"
        SSF_HF = "ssf_HF"
        WVG = "wvg"

    PLOT_X_AXIS = {
        ResultNames.ALPHA.value: ResultNames.FXC_GRID.value,
        ResultNames.ADR.value: ResultNames.WVG.value,
        ResultNames.BF.value: ResultNames.WVG.value,
        ResultNames.IDR.value: ResultNames.WVG.value,
        ResultNames.RDF.value: ResultNames.RDF_GRID.value,
        ResultNames.SDR.value: ResultNames.WVG.value,
        ResultNames.SLFC.value: ResultNames.WVG.value,
        ResultNames.SSF.value: ResultNames.WVG.value,
        ResultNames.SSF_HF.value: ResultNames.WVG.value,
    }

    PLOT_X_LABEL = {
        ResultNames.FXC_GRID.value: "Coupling parameter",
        ResultNames.RDF_GRID.value: "Inter-particle distance",
        ResultNames.WVG.value: "Wave-vector",
    }

    PLOT_Y_LABEL = {
        ResultNames.ALPHA.value: "Free Parameter for VS schemes",
        ResultNames.ADR.value: "Auxiliary density response",
        ResultNames.BF.value: "Bridge function adder",
        ResultNames.IDR.value: "Ideal density response",
        ResultNames.RDF.value: "Radial distribution function",
        ResultNames.SDR.value: "Static density response",
        ResultNames.SLFC.value: "Static local field correction",
        ResultNames.SSF.value: "Static structure factor",
        ResultNames.SSF_HF.value: "Hartree-Fock static structure factor",
    }

    # Read runs in the database
    @staticmethod
    def inspect_runs(database_name: str | None = None) -> dict:
        """Reads runs from the database and returns the content in the form of a dictionary.

        Args:
            database_name (str, optional): Name of the database to read from. Defaults to None.

        Returns:
            A dictionary whose keys are the run ids and values are the corresponding runs information.
        """
        db_handler = database.DataBaseHandler(database_name)
        return db_handler.inspect_runs()

    # Read runs in the database
    @staticmethod
    def read_run(
        run_id: int,
        database_name: str | None = None,
        input_names: list[str] | None = None,
        result_names: list[str] | None = None,
    ) -> dict:
        """Reads runs from the database and returns the content in the form of a dictionary.

        Args:
            database_name (str, optional): Name of the database to read from. Defaults to None.

        Returns:
            A dictionary whose keys are the run ids and values are the corresponding runs information.
        """
        db_handler = database.DataBaseHandler(database_name)
        return db_handler.get_run(run_id, input_names, result_names)

    # Read inputs in the database
    @staticmethod
    def read_inputs(
        run_id: int, database_name: str | None = None, names: list[str] | None = None
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
        return db_handler.get_inputs(run_id, names if names is not None else [])

    # Read results in the database
    @staticmethod
    def read_results(
        run_id: int, database_name: str | None = None, names: list[str] | None = None
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
        return db_handler.get_results(run_id, names)

    # Plot from data in hdf file
    @staticmethod
    def plot(
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
        for y_name in to_plot:
            if y_name not in HDF.PLOT_Y_LABEL.keys():
                raise ValueError(f"Unknown quantity to plot: {y_name}")
            x_name = HDF.PLOT_X_AXIS[y_name]
            if y_name in [HDF.ResultNames.ADR.value, HDF.ResultNames.IDR.value]:
                if matsubara is None:
                    matsubara = HDF.read_inputs(run_id, database_name, "matsubara")
                HDF._plot_1D_parametric(
                    x_name, y_name, matsubara, run_id, database_name
                )
            else:
                HDF._plot_1D(x_name, y_name, run_id, database_name)

    @staticmethod
    def _plot_1D(x_name: str, y_name: str, run_id: id, database_name: str | None):
        data = HDF.read_results(run_id, database_name, [x_name, y_name])
        Plot.plot_1D(
            data[x_name],
            data[y_name],
            HDF.PLOT_X_LABEL.get(x_name, ""),
            HDF.PLOT_Y_LABEL.get(y_name, ""),
        )

    @staticmethod
    def _plot_1D_parametric(
        x_name: str, y_name: str, parameters, run_id: id, database_name: str | None
    ):
        data = HDF.read_results(run_id, database_name, [x_name, y_name])
        Plot.plot_1D_parametric(
            data[x_name],
            data[y_name],
            HDF.PLOT_X_LABEL.get(x_name, ""),
            HDF.PLOT_Y_LABEL.get(y_name, ""),
            parameters,
        )


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
