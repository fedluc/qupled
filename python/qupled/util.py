import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpi4py import MPI as MPINative # This import initializes MPI
import qupled.qupled as qp

# -----------------------------------------------------------------------
# Hdf class
# -----------------------------------------------------------------------

class Hdf():

    """ Class to manipulate the output hdf files produced when a scheme is solved. """
    
    # Construct
    def __init__(self):
        # The first entry is a descriptive name of the 
        self.entries = {
            "adr"        : self.Entries("Auxiliary density response", "numpy2D"),
            "bf"         : self.Entries("Bridge function adder", "numpy"),
            "coupling"   : self.Entries("Coupling parameter", "number"),
            "cutoff"     : self.Entries("Cutoff for the wave-vector grid", "number"),
            "degeneracy" : self.Entries("Degeneracy parameter", "number"),
            "fxcGrid"    : self.Entries("Coupling parameter", "numpy"),
            "fxci"       : self.Entries("Free Energy integrand", "numpy2D"),
            "matsubara"  : self.Entries("Number of matsubara frequencies", "number"),
            "idr"        : self.Entries("Ideal density response", "numpy2D"),
            "resolution" : self.Entries("Resolution for the wave-vector grid", "number"),
            "rdf"        : self.Entries("Radial distribution function", "numpy"),
            "rdfGrid"    : self.Entries("Inter-particle distance", "numpy"), 
            "sdr"        : self.Entries("Static density response", "numpy"),
            "slfc"       : self.Entries("Static local field correction", "numpy"),
            "ssf"        : self.Entries("Static structure factor", "numpy"),
            "ssfHF"      : self.Entries("Hartree-Fock static structure factor", "numpy"),
            "theory"     : self.Entries("Theory that is being solved", "string"),
            "wvg"        : self.Entries("Wave-vector", "numpy"),
        }

    # Structure used to cathegorize the entries stored in the hdf file
    class Entries():
        def __init__(self, description, entryType):
            self.description = description # Descriptive string of the entry
            self.entryType = entryType # Type of entry (numpy, numpy2, number or string)
            assert(self.entryType == "numpy"
                   or self.entryType == "numpy2D"
                   or self.entryType == "number"
                   or self.entryType == "string")
            
    # Read data in hdf file
    def read(self, hdf : str, toRead : list[str]) -> dict:
        """ Reads an hdf file produced by coupled and returns the content in the form of a dictionary

        Args:  
            hdf: Name of the hdf file to read  
            toRead: A list of quantities to read. The list of quantities that can be extracted from
                the hdf file can be obtained by running :func:`~qupled.util.Hdf.inspect`
        
        Returns:
            A dictionary whose entries are the quantities listed in toRead
        
        """
        output = dict.fromkeys(toRead)
        for name in toRead:
            if ( self.entries[name].entryType == "numpy" ) :
                output[name] = pd.read_hdf(hdf, name)[0].to_numpy()
            elif ( self.entries[name].entryType == "numpy2D" ) :
                output[name] = pd.read_hdf(hdf, name).to_numpy()
            elif ( self.entries[name].entryType == "number") :
                output[name] = pd.read_hdf(hdf, "inputs")[name][0].tolist()
            elif ( self.entries[name].entryType == "string") :
                output[name] = pd.read_hdf(hdf, "inputs")[name][0]
            else:
                sys.exit("Unknown entry type")
        return output

    # Get all quantities stored in an hdf file
    def inspect(self, hdf : str) -> dict:
        """ Allows to obtain a summary of the quantities stored in an hdf file

        Args:  
            hdf: Name of the hdf file to inspect
        
        Returns:
            A dictionary containing all the quantities stored in the hdf file and a brief description for
            each quantity
        
        """
        with pd.HDFStore(hdf, mode='r') as store:
            datasetNames =  [name[1:] if name.startswith('/') else name for name in store.keys()]
            if "inputs" in datasetNames:
                datasetNames.remove("inputs")
                for name in store["inputs"].keys() : datasetNames.append(name)
        output = dict.fromkeys(datasetNames)
        for key in output.keys():
            output[key] = self.entries[key].description;
        return output
        
    # Plot from data in hdf file
    def plot(self, hdf : str, toPlot : list[str], matsubara : np.array = None) -> None:
        """ Plots the results stored in an hdf file.

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
        for name in toPlot :
            if ( name == "rdf" ) :
                x = self.read(hdf, [name, "rdfGrid"])
                Plot.plot1D(x["rdfGrid"], x[name], self.entries["rdfGrid"].description, self.entries[name].description)
            elif ( name in ["adr", "idr"] ) :
                x = self.read(hdf, [name, "wvg", "matsubara"])
                if (matsubara is None) : matsubara = np.arange(x["matsubara"])
                Plot.plot1DParametric(x["wvg"], x[name], self.entries["wvg"].description, self.entries[name].description, matsubara)
            elif ( name == "fxci" ) :
                x = self.read(hdf, [name, "fxcGrid"])
                Plot.plot1D(x["fxcGrid"], x[name][1,:], self.entries["fxcGrid"].description, self.entries[name].description)
            elif ( name in ["sdr", "slfc", "ssf", "ssfHF"] ) :
                x = self.read(hdf, [name, "wvg"])
                Plot.plot1D(x["wvg"], x[name], self.entries["wvg"].description, self.entries[name].description)
            else:
                sys.exit("Unknown quantity to plot")

        
    def computeRdf(self, hdf, rdfGrid : np.array = None, saveRdf : bool = True) -> None:
        """ Computes the radial distribution function and returns it as a numpy array.

        Args:  
            hdf: Name of the hdf file to load the structural properties from
            rdfGrid: A numpy array specifing the grid used to compute the radial distribution function
                (default = None, i.e. rdfGrid = np.arange(0.01, 10.0, 0.01))
            saveRdf: Flag marking whether the rdf data should be added to the hdf file (default = True)

        Returns:
            The radial distribution function
        
        """
        hdfData = self.read(hdf, ["wvg", "ssf"])
        if (rdfGrid is None) : rdfGrid = np.arange(0.01, 10.0, 0.01)
        rdf = qp.computeRdf(rdfGrid, hdfData["wvg"], hdfData["ssf"])
        if (saveRdf):
            pd.DataFrame(rdfGrid).to_hdf(hdf, key="rdfGrid", mode="r+")
            pd.DataFrame(rdf).to_hdf(hdf, key="rdf", mode="r+")
        return rdf

    def computeInternalEnergy(self, hdf : str) -> float:
        """ Computes the internal energy and returns it to output.

        Args:  
            hdf: Name of the hdf file to load the structural properties from

        Returns:
            The internal energy
        """
        hdfData = self.read(hdf, ["wvg", "ssf", "coupling"])
        return qp.computeInternalEnergy(hdfData["wvg"], hdfData["ssf"], hdfData["coupling"])
    
# -----------------------------------------------------------------------
# Plot class
# -----------------------------------------------------------------------
    
class Plot():

    """ Class to collect methods used for plotting """
    
    # One dimensional plots
    def plot1D(x, y, xlabel, ylabel):
        """ Produces the plot of a one-dimensional quantity.

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
        """ Produces the plot of a one-dimensional quantity that depends on an external parameter.
        
        Positional arguments:  
        x -- data for the x-axis (a numpy array)  
        y -- data for the y-axis (a two-dimensional numpy array)  
        xlabel -- label for the x-axis (a string)  
        ylabel -- label for the y-axis (a string)  
        parameters -- list of parameters for which the results should be plotted
        """
        numParameters = parameters.size
        cmap = cm.get_cmap(name="viridis")
        for i in np.arange(numParameters):
            color = cmap(1.0*i/numParameters)
            plt.plot(x, y[:,parameters[i]], color=color)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()        


# -----------------------------------------------------------------------
# MPI class
# -----------------------------------------------------------------------

class MPI():

    """ Class to handle the calls to the MPI API """
    def __init__(self):
        self.communicator = MPINative.COMM_WORLD


    def getRank(self):
        return self.communicator.Get_rank()

    def isRoot(self):
        return self.getRank() == 0

    def barrier(self):
        self.communicator.Barrier()
