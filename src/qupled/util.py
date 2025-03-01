from enum import Enum
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
import functools
from qupled import native

# -----------------------------------------------------------------------
# Hdf class
# -----------------------------------------------------------------------


class Hdf:
    """Class to manipulate the output hdf files produced when a scheme is solved."""

    class EntryKeys(Enum):
        ALPHA = "alpha"
        ADR = "adr"
        BF = "bf"
        COUPLING = "coupling"
        CUTOFF = "cutoff"
        FREQUENCY_CUTOFF = "frequencyCutoff"
        DEGENERACY = "degeneracy"
        ERROR = "error"
        FXC_GRID = "fxcGrid"
        FXCI = "fxci"
        MATSUBARA = "matsubara"
        IDR = "idr"
        INFO = "info"
        RESOLUTION = "resolution"
        RDF = "rdf"
        RDF_GRID = "rdfGrid"
        SDR = "sdr"
        SLFC = "slfc"
        SSF = "ssf"
        SSF_HF = "ssfHF"
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
            Hdf.EntryKeys.ALPHA.value: self.Entries(
                "Free Parameter for VS schemes", Hdf.EntryType.NUMPY.value
            ),
            Hdf.EntryKeys.ADR.value: self.Entries(
                "Auxiliary density response", Hdf.EntryType.NUMPY2D.value
            ),
            Hdf.EntryKeys.BF.value: self.Entries(
                "Bridge function adder", Hdf.EntryType.NUMPY.value
            ),
            Hdf.EntryKeys.COUPLING.value: self.Entries(
                "Coupling parameter", Hdf.EntryType.NUMBER.value
            ),
            Hdf.EntryKeys.CUTOFF.value: self.Entries(
                "Cutoff for the wave-vector grid", Hdf.EntryType.NUMBER.value
            ),
            Hdf.EntryKeys.FREQUENCY_CUTOFF.value: self.Entries(
                "Cutoff for the frequency", Hdf.EntryType.NUMBER.value
            ),
            Hdf.EntryKeys.DEGENERACY.value: self.Entries(
                "Degeneracy parameter", Hdf.EntryType.NUMBER.value
            ),
            Hdf.EntryKeys.ERROR.value: self.Entries(
                "Residual error in the solution", Hdf.EntryType.NUMBER.value
            ),
            Hdf.EntryKeys.FXC_GRID.value: self.Entries(
                "Coupling parameter", Hdf.EntryType.NUMPY.value
            ),
            Hdf.EntryKeys.FXCI.value: self.Entries(
                "Free Energy integrand", Hdf.EntryType.NUMPY2D.value
            ),
            Hdf.EntryKeys.MATSUBARA.value: self.Entries(
                "Number of matsubara frequencies", Hdf.EntryType.NUMBER.value
            ),
            Hdf.EntryKeys.IDR.value: self.Entries(
                "Ideal density response", Hdf.EntryType.NUMPY2D.value
            ),
            Hdf.EntryKeys.RESOLUTION.value: self.Entries(
                "Resolution for the wave-vector grid", Hdf.EntryType.NUMBER.value
            ),
            Hdf.EntryKeys.RDF.value: self.Entries(
                "Radial distribution function", Hdf.EntryType.NUMPY.value
            ),
            Hdf.EntryKeys.RDF_GRID.value: self.Entries(
                "Inter-particle distance", Hdf.EntryType.NUMPY.value
            ),
            Hdf.EntryKeys.SDR.value: self.Entries(
                "Static density response", Hdf.EntryType.NUMPY.value
            ),
            Hdf.EntryKeys.SLFC.value: self.Entries(
                "Static local field correction", Hdf.EntryType.NUMPY.value
            ),
            Hdf.EntryKeys.SSF.value: self.Entries(
                "Static structure factor", Hdf.EntryType.NUMPY.value
            ),
            Hdf.EntryKeys.SSF_HF.value: self.Entries(
                "Hartree-Fock static structure factor", Hdf.EntryType.NUMPY.value
            ),
            Hdf.EntryKeys.THEORY.value: self.Entries(
                "Theory that is being solved", Hdf.EntryType.STRING.value
            ),
            Hdf.EntryKeys.WVG.value: self.Entries(
                "Wave-vector", Hdf.EntryType.NUMPY.value
            ),
        }

    # Structure used to cathegorize the entries stored in the hdf file
    class Entries:
        def __init__(self, description, entryType):
            self.description = description  # Descriptive string of the entry
            self.entryType = (
                entryType  # Type of entry (numpy, numpy2, number or string)
            )
            assert (
                self.entryType == Hdf.EntryType.NUMPY.value
                or self.entryType == Hdf.EntryType.NUMPY2D.value
                or self.entryType == Hdf.EntryType.NUMBER.value
                or self.entryType == Hdf.EntryType.STRING.value
            )

    # Read data in hdf file
    def read(self, hdf: str, toRead: list[str]) -> dict:
        """Reads an hdf file produced by coupled and returns the content in the form of a dictionary

        Args:
            hdf: Name of the hdf file to read
            toRead: A list of quantities to read. The list of quantities that can be extracted from
                the hdf file can be obtained by running :func:`~qupled.util.Hdf.inspect`

        Returns:
            A dictionary whose entries are the quantities listed in toRead

        """
        output = dict.fromkeys(toRead)
        with pd.HDFStore(hdf, mode="r") as store:
            for name in toRead:
                if name not in self.entries:
                    raise KeyError(f"Unknown entry: {name}")
                if self.entries[name].entryType == Hdf.EntryType.NUMPY.value:
                    output[name] = store[name][0].to_numpy()
                elif self.entries[name].entryType == Hdf.EntryType.NUMPY2D.value:
                    output[name] = store[name].to_numpy()
                elif self.entries[name].entryType == Hdf.EntryType.NUMBER.value:
                    output[name] = (
                        store[Hdf.EntryKeys.INFO.value][name].iloc[0].tolist()
                    )
                elif self.entries[name].entryType == Hdf.EntryType.STRING.value:
                    output[name] = store[Hdf.EntryKeys.INFO.value][name].iloc[0]
                else:
                    raise ValueError(
                        f"Unknown entry type: {self.entries[name].entryType}"
                    )
        return output

    # Get all quantities stored in an hdf file
    def inspect(self, hdf: str) -> dict:
        """Allows to obtain a summary of the quantities stored in an hdf file

        Args:
            hdf: Name of the hdf file to inspect

        Returns:
            A dictionary containing all the quantities stored in the hdf file and a brief description for
            each quantity

        """
        with pd.HDFStore(hdf, mode="r") as store:
            datasetNames = [
                name[1:] if name.startswith("/") else name for name in store.keys()
            ]
            if Hdf.EntryKeys.INFO.value in datasetNames:
                datasetNames.remove(Hdf.EntryKeys.INFO.value)
                for name in store[Hdf.EntryKeys.INFO.value].keys():
                    datasetNames.append(name)
        output = dict.fromkeys(datasetNames)
        for key in output.keys():
            output[key] = self.entries[key].description
        return output

    # Plot from data in hdf file
    def plot(self, hdf: str, toPlot: list[str], matsubara: np.array = None) -> None:
        """Plots the results stored in an hdf file.

        Args:
            hdf: Name of the hdf file
            toPlot: A list of quantities to plot. Allowed quantities include adr (auxiliary density response),
                bf (bridge function adder), fxci (free energy integrand), idr (ideal density response), rdf
                (radial distribution function), sdr (static density response), slfc (static local field correction)
                ssf (static structure factor) and ssfHF (Hartree-Fock static structure factor).
                If the hdf file does not contain the specified quantity, an error is thrown
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Defaults to  None, all matsubara frequencies are plotted)

        """
        for name in toPlot:
            description = (
                self.entries[name].description if name in self.entries.keys() else ""
            )
            if name == Hdf.EntryKeys.RDF.value:
                x = self.read(hdf, [name, Hdf.EntryKeys.RDF_GRID.value])
                Plot.plot1D(
                    x[Hdf.EntryKeys.RDF_GRID.value],
                    x[name],
                    self.entries[Hdf.EntryKeys.RDF_GRID.value].description,
                    description,
                )
            elif name in [Hdf.EntryKeys.ADR.value, Hdf.EntryKeys.IDR.value]:
                x = self.read(
                    hdf, [name, Hdf.EntryKeys.WVG.value, Hdf.EntryKeys.MATSUBARA.value]
                )
                if matsubara is None:
                    matsubara = np.arange(x[Hdf.EntryKeys.MATSUBARA.value])
                Plot.plot1DParametric(
                    x[Hdf.EntryKeys.WVG.value],
                    x[name],
                    self.entries[Hdf.EntryKeys.WVG.value].description,
                    description,
                    matsubara,
                )
            elif name == Hdf.EntryKeys.FXCI.value:
                x = self.read(hdf, [name, Hdf.EntryKeys.FXC_GRID.value])
                Plot.plot1D(
                    x[Hdf.EntryKeys.FXC_GRID.value],
                    x[name][1, :],
                    self.entries[Hdf.EntryKeys.FXC_GRID.value].description,
                    description,
                )
            elif name in [
                Hdf.EntryKeys.BF.value,
                Hdf.EntryKeys.SDR.value,
                Hdf.EntryKeys.SLFC.value,
                Hdf.EntryKeys.SSF.value,
                Hdf.EntryKeys.SSF_HF.value,
            ]:
                x = self.read(hdf, [name, Hdf.EntryKeys.WVG.value])
                Plot.plot1D(
                    x[Hdf.EntryKeys.WVG.value],
                    x[name],
                    self.entries[Hdf.EntryKeys.WVG.value].description,
                    self.entries[name].description,
                )
            elif name == Hdf.EntryKeys.ALPHA.value:
                x = self.read(hdf, [name, Hdf.EntryKeys.FXC_GRID.value])
                Plot.plot1D(
                    x[Hdf.EntryKeys.FXC_GRID.value][::2],
                    x[name][::2],
                    self.entries[Hdf.EntryKeys.FXC_GRID.value].description,
                    self.entries[name].description,
                )
            else:
                raise ValueError(f"Unknown quantity to plot: {name}")

    def computeRdf(
        self, hdf: str, rdfGrid: np.array = None, saveRdf: bool = True
    ) -> None:
        """Computes the radial distribution function and returns it as a numpy array.

        Args:
            hdf: Name of the hdf file to load the structural properties from
            rdfGrid: A numpy array specifing the grid used to compute the radial distribution function
                (default = None, i.e. rdfGrid = np.arange(0.0, 10.0, 0.01))
            saveRdf: Flag marking whether the rdf data should be added to the hdf file (default = True)

        Returns:
            The radial distribution function

        """
        hdfData = self.read(hdf, [Hdf.EntryKeys.WVG.value, Hdf.EntryKeys.SSF.value])
        if rdfGrid is None:
            rdfGrid = np.arange(0.0, 10.0, 0.01)
        rdf = native.compute_rdf(
            rdfGrid, hdfData[Hdf.EntryKeys.WVG.value], hdfData[Hdf.EntryKeys.SSF.value]
        )
        if saveRdf:
            pd.DataFrame(rdfGrid).to_hdf(
                hdf, key=Hdf.EntryKeys.RDF_GRID.value, mode="r+"
            )
            pd.DataFrame(rdf).to_hdf(hdf, key=Hdf.EntryKeys.RDF.value, mode="r+")
        return rdf

    def computeInternalEnergy(self, hdf: str) -> float:
        """Computes the internal energy and returns it to output.

        Args:
            hdf: Name of the hdf file to load the structural properties from

        Returns:
            The internal energy
        """
        hdfData = self.read(
            hdf,
            [
                Hdf.EntryKeys.WVG.value,
                Hdf.EntryKeys.SSF.value,
                Hdf.EntryKeys.COUPLING.value,
            ],
        )
        return native.compute_internal_energy(
            hdfData[Hdf.EntryKeys.WVG.value],
            hdfData[Hdf.EntryKeys.SSF.value],
            hdfData[Hdf.EntryKeys.COUPLING.value],
        )


# -----------------------------------------------------------------------
# Plot class
# -----------------------------------------------------------------------


class Plot:
    """Class to collect methods used for plotting"""

    # One dimensional plots
    def plot1D(x, y, xlabel, ylabel):
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

    # One dimensional plots with one parameter"
    def plot1DParametric(x, y, xlabel, ylabel, parameters):
        """Produces the plot of a one-dimensional quantity that depends on an external parameter.

        Positional arguments:
        x -- data for the x-axis (a numpy array)
        y -- data for the y-axis (a two-dimensional numpy array)
        xlabel -- label for the x-axis (a string)
        ylabel -- label for the y-axis (a string)
        parameters -- list of parameters for which the results should be plotted
        """
        numParameters = parameters.size
        cmap = cm["viridis"]
        for i in np.arange(numParameters):
            color = cmap(1.0 * i / numParameters)
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
        self.qpMPI = native.MPI

    def getRank(self):
        """Get rank of the process"""
        return self.qpMPI.rank()

    def isRoot(self):
        """Check if the current process is root (rank 0)"""
        return self.qpMPI.is_root()

    def barrier(self):
        """Setup and MPI barrier"""
        self.qpMPI.barrier()

    def timer(self):
        """Get wall time"""
        return self.qpMPI.timer()

    @staticmethod
    def runOnlyOnRoot(func):
        """Python decorator for all methods that have to be run only by root"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if MPI().isRoot():
                return func(*args, **kwargs)

        return wrapper

    @staticmethod
    def synchronizeRanks(func):
        """Python decorator for all methods that need to rank synchronization"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            func(*args, **kwargs)
            MPI().barrier()

        return wrapper

    @staticmethod
    def recordTime(func):
        """Python decorator for all methods that have to be timed"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            tic = MPI().timer()
            func(*args, **kwargs)
            toc = MPI().timer()
            dt = toc - tic
            hours = dt // 3600
            minutes = (dt % 3600) // 60
            seconds = dt % 60
            if MPI().isRoot():
                if hours > 0:
                    print("Elapsed time: %d h, %d m, %d s." % (hours, minutes, seconds))
                elif minutes > 0:
                    print("Elapsed time: %d m, %d s." % (minutes, seconds))
                else:
                    print("Elapsed time: %.1f s." % seconds)

        return wrapper
