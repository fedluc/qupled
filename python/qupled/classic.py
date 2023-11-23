import sys
import os
from shutil import rmtree
from glob import glob
import zipfile as zf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import qupled.qupled as qp

class Stls():

    """Class to solve the STLS scheme.

    Class used to setup and solve the classical STLS scheme as described by
    `Tanaka and Ichimaru <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_.
    The inputs used to solve the scheme are defined when creating the class, but can be
    later modified by changing the attribute :obj:`inputs`. After the solution is completed
    the results are saved to an hdf file and can be plotted via the method :obj:`plot`

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.  
        chemicalPotential: Initial guess for the chemical potential, defaults to [-10.0, 10.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for convergence, defaults to 1.0e-5.
        mixing: Mixing parameter for iterative solution, defaults to 1.0.  
        guess:  Initial guess for the iterative solution, defaults to None, i.e. slfc = 0.
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
    """
    
    # Constructor
    def __init__(self,
                 coupling : float,
                 degeneracy : float,
                 chemicalPotential : list[float] = [-10.0,10.0],
                 cutoff : float = 10.0,
                 error : float = 1.0e-5,
                 mixing : float = 1.0,
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1 ):
        # Allowed theories
        self.allowedTheories = ["STLS"]
        # Input object
        self.inputs : qupled.qupled.StlsInput = qp.StlsInput() #: Inputs to solve the scheme.
        self._setInputs(coupling, degeneracy, "STLS", chemicalPotential,
                        cutoff, error, mixing, iterations, matsubara,
                        outputFrequency, recoveryFile, resolution)
        # Scheme to solve and associated input and solution
        self.scheme : qp.Stls = None
        self.solutionAvailable = False
        # File to store output on disk
        self.hdfFileName : str = None

    # Setup inputs object
    def _setInputs(self,
                   coupling : float,
                   degeneracy : float,
                   theory : str,
                   chemicalPotential : list[float],
                   cutoff : float,
                   error : float,
                   mixing : float,
                   iterations : int,
                   matsubara : int,
                   outputFrequency : int,
                   recoveryFile : str,
                   resolution : float) -> None:
        """ Sets up the content of :obj:`inputs` """
        self.inputs.coupling = coupling
        self.inputs.degeneracy = degeneracy
        self.inputs.theory = theory
        self.inputs.chemicalPotential = chemicalPotential
        self.inputs.cutoff = cutoff
        self.inputs.error = error
        # if (guess is not None): self.inputs.guess = guess
        self.inputs.mixing = mixing
        self.inputs.iterations = iterations
        self.inputs.matsubara = matsubara
        self.inputs.outputFrequency = outputFrequency
        if (recoveryFile is not None): self.inputs.recoveryFile = recoveryFile
        self.inputs.resolution = resolution
        self.inputs.intError = 1.0e-5
        self.inputs.threads = 1
        
    # Check input before computing
    def _checkInputs(self) -> None:
        """ Checks that the content of :obj:`inputs` is correct """
        if (self.inputs.theory not in self.allowedTheories):
            sys.exit("Invalid dielectric theory")

    # Compute
    def compute(self) -> None:
        """ Solves the scheme and saves the results to and hdf file. See the method :func:`~qupled.classic.Stls.save`
        to see which results are saved
        """
        self._checkInputs()
        self.scheme = qp.Stls(self.inputs)
        status = self.scheme.compute()
        self._checkStatusAndClean(status)        
        self._setHdfFile()
        self._save()
        self.solutionAvailable = True
        
    # Check that the dielectric scheme was solved without errors
    def _checkStatusAndClean(self, status : bool) -> None:
        """ Checks that the scheme was solved correctly and removes temporarary files generated at run-time
        
           Args:
               status: status obtained from the native code. If status == 0 the scheme was solved correctly.
        """
        if (status == 0):
            if os.path.isfile(self.scheme.recovery) : os.remove(self.scheme.recovery)
            print("Dielectric theory solved successfully!")
        else:
            sys.exit("Error while solving the dielectric theory")
    
    # Save results to disk
    def _setHdfFile(self) -> None:
        """ Sets the name of the hdf file used to store the output """
        self.hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (self.inputs.coupling,
                                                         self.inputs.degeneracy,
                                                         self.inputs.theory)
    
    def _save(self) -> None:
        """ Stores the results obtained by solving the scheme.

        The results are stored as pandas dataframes in an hdf file with the following keywords:
        
        - inputs: A dataframe containing information on the input parameters, it includes:
        
          - coupling: the coupling parameter,
          - degeneracy: the degeneracy parameter,
          - theory: the theory that is being solved,  
          - error: the minimum error for convergence,   
          - resolution: the resolution in the wave-vector grid, 
          - cutoff: the cutoff in the wave-vector grid, 
          - matsubara: the number of matsubara frequencies
        
        - idr (*ndarray*, 2D): the ideal density response  
        - sdr (*ndarray*):  the static density response   
        - slfc (*ndarray*):  the static local field correction   
        - ssf (*ndarray*):  the static structure factor   
        - ssfHF (*ndarray*):  the Hartree-Fock static structure factor 
        - wvg (*ndarray*):  the wave-vector grid 
        
        If the radial distribution function was computed (see computeRdf), then the hdf file contains
        two additional keywords:
        
        - rdf (*ndarray*):  the radial distribution function   
        - rdfGrid (*ndarray*):  the grid used to compute the radial distribution function 
        
        """
        pd.DataFrame({
            "coupling" : self.inputs.coupling,
            "degeneracy" : self.inputs.degeneracy,
            "theory" : self.inputs.theory,
            "error" : self.inputs.error,
            "resolution" : self.inputs.resolution,
            "cutoff" : self.inputs.cutoff,
            "matsubara" : self.inputs.matsubara
            }, index=["inputs"]).to_hdf(self.hdfFileName, key="inputs", mode="w")
        pd.DataFrame(self.scheme.idr).to_hdf(self.hdfFileName, key="idr")
        pd.DataFrame(self.scheme.sdr).to_hdf(self.hdfFileName, key="sdr")
        pd.DataFrame(self.scheme.slfc).to_hdf(self.hdfFileName, key="slfc")
        pd.DataFrame(self.scheme.ssf).to_hdf(self.hdfFileName, key="ssf")
        pd.DataFrame(self.scheme.ssfHF).to_hdf(self.hdfFileName, key="ssfHF")
        pd.DataFrame(self.scheme.wvg).to_hdf(self.hdfFileName, key="wvg")
        
    # Compute radial distribution function
    def computeRdf(self, rdfGrid : np.ndarray, writeToHdf : bool = True) -> np.array:
        """ Computes the radial distribution function from the data stored in :obj:`scheme`.
        
        Args:
            rdfGrid: The grid used to compute the radial distribution function
            writeToHdf: Flag marking whether the rdf data should be added to the output hdf file, default to True

        Returns:
            The radial distribution function
        
        """
        if (not self.solutionAvailable):
            sys.exit("No solution to compute the radial distribution function")
        rdf = self.scheme.rdf(rdfGrid)
        if (writeToHdf) :
            pd.DataFrame(rdfGrid).to_hdf(self.hdfFileName, key="rdfGrid", mode="r+")
            pd.DataFrame(rdf).to_hdf(self.hdfFileName, key="rdf", mode="r+")
        return rdf

    # Compute the internal energy
    def computeInternalEnergy(self) -> float:
        """ Computes the internal energy from the data stored in :obj:`scheme`.

        Returns:
            The internal energy
        
        """
        if (not self.solutionAvailables):
            sys.exit("No solution to compute the internal energy")
        return qp.computeInternalEnergy(self.scheme.wvg, self.scheme.ssf, self.inputs.coupling)
        
    # Plot results
    def plot(self, toPlot : list[str], matsubara : np.ndarray = None, rdfGrid : np.ndarray = None) -> None:
        """ Plots the results obtained stored in :obj:`scheme`.

        Args:  
            toPlot: A list of quantities to plot. This list can include idr (ideal density response), rdf
                (radial distribution function), sdr (static density response), slfc (static local field correction)
                ssf (static structure factor) and ssfHF (Hartree-Fock static structure factor)  
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Default = None, see :func:`~qupled.classic.Stls.plotIdr`)  
            rdfGrid: The grid used to compute the radial distribution function. Applies only when the radial
                distribution function is plotted (Default = None, see :func:`~qupled.classic.Stls.computeRdf`)
        
        """
        wvg = self.scheme.wvg
        xlabel = "Wave vector"
        if ("idr" in toPlot):
            self._plotIdr(matsubara)
        if ("rdf" in toPlot):
            self._plotRdf(rdfGrid)
        if ("sdr" in toPlot):
            Plot.plot1D(wvg, self.scheme.sdr, xlabel, "Static density response")
        if ("slfc" in toPlot):
            Plot.plot1D(wvg, self.scheme.slfc, xlabel, "Static local field correction")
        if ("ssf" in toPlot):
            Plot.plot1D(wvg, self.scheme.ssf, xlabel, "Static structure factor")
        if ("ssfHF" in toPlot):
            Plot.plot1D(wvg, self.scheme.ssfHF, xlabel, "Hartree-Fock static structure factor")
        
    def _plotIdr(self, matsubara : np.ndarray = None) -> None:
        """ Plots the ideal density response.
        
        Args:  
            matsubara:  A list of matsubara frequencies to plot. (Default =  all matsubara frequencies are plotted)
        
        """
        if (self.inputs.degeneracy == 0) : return
        if (matsubara is None) : matsubara = np.arange(self.inputs.matsubara)
        Plot.plot1DParametric(self.scheme.wvg, self.scheme.idr,
                              "Wave vector", "Ideal density response",
                              matsubara)
        
    def _plotRdf(self, rdfGrid : np.ndarray = None) -> None:
        """ Plot the radial distribution function.
        
        Args:  
            rdfGrid: The grid used to compute the radial distribution function. Applies only when the radial
                distribution function is plotted (Default = None, i.e.  numpy.arange(0.01, 10.0, 0.01)`)
          
        """
        if (rdfGrid is None) : rdfGrid = np.arange(0.01, 10.0, 0.01)
        rdf = self.computeRdf(rdfGrid)
        Plot.plot1D(rdfGrid, rdf, "Inter-particle distance", "radial distribution function")
        
        
class StlsIet(Stls):

    """Class to solve the STLS-IET schemes.

    Class used to setup and solve the classical STLS-IET scheme as described by
    `Tanaka <https://aip.scitation.org/doi/full/10.1063/1.4969071>`_ and by
    `Tolias and collaborators <https://aip.scitation.org/doi/full/10.1063/1.4969071>`_.
    This class inherits most of its methods and attributes from :obj:`~qupled.classic.Stls`

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.  
        chemicalPotential: Initial guess for the chemical potential, defaults to [-10.0, 10.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for convergence, defaults to 1.0e-5.
        mapping: Classical to quantum mapping. See :func:`~qupled.qupled.StlsInput.iet`
        mixing: Mixing parameter for iterative solution, defaults to 1.0.  
        guess:  Initial guess for the iterative solution, defaults to None, i.e. slfc = 0.
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
        scheme2DIntegrals: numerical scheme used to solve two-dimensional integrals. See :func:`~qupled.qupled.Input.int2DScheme`
    """
    
    # Constructor
    def __init__(self,
                 coupling : float,
                 degeneracy : float,
                 theory : str,
                 chemicalPotential : list[float] = [-10.0,10.0],
                 cutoff : float = 10.0,
                 error : float = 1.0e-5,
                 mapping : str = "standard",
                 mixing : float = 1.0,
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1,
                 scheme2DIntegrals : str = "full"):
        # Allowed theories
        self.allowedTheories : ["STLS-HNC", "STLS-IOI", "STLS-LCT"]
        # Input object
        self.inputs : qupled.qupled.StlsInput = qp.StlsInput() #: Inputs to solve the scheme.
        super()._setInputs(coupling, degeneracy, theory, chemicalPotential,
                           cutoff, error, mixing, iterations, matsubara,
                           outputFrequency, recoveryFile, resolution)
        self.inputs.iet = mapping
        self.inputs.int2DScheme = scheme2DIntegrals
        self._checkInputs()
        # Scheme to solve and associated input and solution
        self.scheme : qp.Stls = None
        # File to store output on disk
        self.hdfFileName = None
            
    # Plot results
    def plot(self, toPlot, matsubara : list[int] = None, rdfGrid : np.ndarray= None) -> None:
        """ Plots the results obtained stored in :obj:`~qupled.classic.Stls.scheme`. Extends 
        :func:`~qupled.classic.Stls.plot` by adding the option to plot the bridge function
        adder by passing `bf` to toPlot
        """
        super().plot(toPlot, matsubara, rdfGrid)
        if ("bf" in toPlot):
            Plot.plot1D(self.scheme.wvg, self.scheme.bf, "Wave vector", "Bridge function adder")
        
    # Save results to disk
    def _save(self) -> None:
        """ Stores the results obtained by solving the scheme. Extends :func:`~qupled.classic.Stls.save`
        by adding the option to save the bridge function adder as a new dataframe in the hdf file. The
        bridge function adder dataframe can be accessed as `bf`
        """
        super()._save()
        pd.DataFrame(self.scheme.bf).to_hdf(self.hdfFileName, key="bf")


class VSStls(Stls):

    """Class to solve the VS-STLS scheme.

    Class used to setup and solve the classical VS-STLS scheme as described by
    `Vashishta and Singwi <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.6.875>`_ and by
    `Sjostrom and Dufty <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.115123>`_.
    This class inherits most of its methods and attributes from :obj:`~qupled.classic.Stls`

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.  
        chemicalPotential: Initial guess for the chemical potential, defaults to [-100.0, 100.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for convergence, defaults to 1.0e-5.
        mixing: Mixing parameter for iterative solution, defaults to 1.0.  
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
        alpha: Initial guess for the free parameter, defaults to [0.5, 1.0]
        couplingResolution: Resolution of the coupling parameter grid, defaults to 0.01
        degeneracyResolution: Resolution of the degeneracy parameter grid, defaults to 0.01
        errorAlpha: Minimum error for convergence in the free parameter iterations, defaults to 1.0e-3
        iterationsAlpha: Maximum number of iterations for the free parameter, defaults to 50
        errorIntegrals: Accuracy (as a relative error) for the integral computations, defaults to 1.0-5
        threads: number of OMP threads for parallel calculations, defualts to 1
    """
    
    # Constructor
    def __init__(self,
                 coupling : float,
                 degeneracy : float,
                 chemicalPotential : list[float] = [-100.0,100.0],
                 cutoff : float = 10.0,
                 error : float = 1.0e-5,
                 mixing : float = 1.0,
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1,
                 alpha : list[float] = [0.5, 1.0],
                 couplingResolution : float = 0.01,
                 degeneracyResolution : float = 0.01,
                 errorAlpha : float = 1.0e-3,
                 iterationsAlpha : int = 50,
                 errorIntegrals : float = 1.0e-5,
                 threads : int = 1):
        # Allowed theories
        self.allowedTheories : list[str] = ["VSSTLS"]
        # Input object
        self.inputs : qupled.qupled.VSStlsInput = qp.VSStlsInput()  #: Inputs to solve the scheme.
        super()._setInputs(coupling, degeneracy, "VSSTLS", chemicalPotential,
                           cutoff, error, mixing, iterations, matsubara,
                           outputFrequency, recoveryFile, resolution)
        self.inputs.alpha = alpha
        self.inputs.couplingResolution = couplingResolution
        self.inputs.degeneracyResolution = degeneracyResolution
        self.inputs.errorAlpha = errorAlpha
        self.inputs.iterationsAlpha = iterationsAlpha
        self.inputs.threads = threads
        self.inputs.intError = errorIntegrals
        # Scheme to solve and associated input and solution
        self.scheme : qp.VSStls = None
        # File to store output on disk
        self.hdfFileName = None

        
    # Compute
    def compute(self) -> None:
        """ Solves the scheme and saves the results to and hdf file. See the method :func:`~qupled.classic.VSStls.save`
        to see which results are saved
        """
        self._checkInputs()
        self.scheme = qp.VSStls(self.inputs)
        status = self.scheme.compute()
        self._checkStatusAndClean(status)        
        self._setHdfFile()
        self._save()

    # Save results
    def _save(self) -> None:
        """ Stores the results obtained by solving the scheme. Extends :func:`~qupled.classic.Stls.save`
        by adding the option to save the free energy integrand and the corresponding coupling parameter grid
        as a new dataframe in the hdf file. The free energy integrand dataframe can be accessed as `fxci`
        and the corresponding coupling parameter grid data frame as `fxcGrid`
        """
        super()._save()
        pd.DataFrame(self.scheme.freeEnergyGrid).to_hdf(self.hdfFileName, key="fxcGrid")
        pd.DataFrame(self.scheme.freeEnergyIntegrand).to_hdf(self.hdfFileName, key="fxci")

    # Plot results        
    def plot(self, toPlot, matsubara : list[int] = None, rdfGrid : np.ndarray= None) -> None:
        """ Plots the results obtained stored in :obj:`~qupled.classic.VSStls.scheme`. Extends 
        :func:`~qupled.classic.Stls.plot` by adding the option to plot the free energy
        integrand by passing `fxci` to toPlot
        """
        super().plot(toPlot, matsubara, rdfGrid)
        if ("fxci" in toPlot):
            Plot.plot1D(self.scheme.freeEnergyGrid, self.scheme.freeEnergyIntegrand[1,:], "Coupling parameter", "Free energy integrand")
        
class Hdf():

    """ Class to manipulate the hdf files produced when a scheme is solved """
            
    # Plot from data in hdf file
    def plot(hdf, toPlot, matsubara = None, rdfGrid = None):
        """ Plots the results stored in an hdf file.

        Positional arguments:  
        hdf -- Name of the hdf file
        toPlot -- A list of quantities to plot. This list can include adr (auxiliary density response),
        bf (bridge function adder), idr (ideal density response), rdf
        (radial distribution function), sdr (static density response), slfc (static local field correction)
        ssf (static structure factor) and ssfHF (Hartree-Fock static structure factor). If the hdf file does not
        contain the specified quantity, an error is thrown
        
        Keyword arguments:  
        matsubara -- A list of matsubara frequencies to plot. Applies only when the idr is plotted.  
        (Default = results for all matsubara frequencies are plotted)    
        rdfGrid -- A numpy array specifing the grid used to compute the radial distribution function
        (Default = None, see computeRdf)  
        
        """
        if (matsubara is None) : matsubara = np.arange(pd.read_hdf(hdf, "inputs")["matsubara"][0].tolist())
        wvg = pd.read_hdf(hdf, "wvg")[0].to_numpy()
        xlabel = "Wave vector"
        if ("adr" in toPlot):
            adr = pd.read_hdf(hdf, "adr").to_numpy()
            Plot.plot1DParametric(wvg, adr, xlabel, "Auxiliary density response", matsubara)
        if ("bf" in toPlot):
            bf = pd.read_hdf(hdf, "bf")[0].to_numpy()
            Plot.plot1D(wvg, bf, xlabel, "Bridge function adder")
        if ("idr" in toPlot):
            idr = pd.read_hdf(hdf, "idr").to_numpy()
            Plot.plot1DParametric(wvg, idr, xlabel, "Ideal density response", matsubara)
        if ("rdf" in toPlot):
            rdf = pd.read_hdf(hdf, "rdf")[0].to_numpy()
            spg = pd.read_hdf(hdf, "rdfGrid")[0].to_numpy()
            Plot.plot1D(spg, rdf, "Inter-particle distance", "Radial distribution function")
        if ("sdr" in toPlot):
            sdr = pd.read_hdf(hdf, "sdr")[0].to_numpy()
            Plot.plot1D(wvg, sdr, xlabel, "Static density response")
        if ("slfc" in toPlot):
            slfc = pd.read_hdf(hdf, "slfc")[0].to_numpy()
            Plot.plot1D(wvg, slfc, xlabel, "Static density response")
        if ("ssf" in toPlot):
            ssf = pd.read_hdf(hdf, "ssf")[0].to_numpy()
            Plot.plot1D(wvg, ssf, xlabel, "Static structure factor")
        if ("ssfHF" in toPlot):
            ssfHF = pd.read_hdf(hdf, "ssfHF")[0].to_numpy()
            Plot.plot1D(wvg, ssfHF, xlabel, "Hartree-Fock static structure factor")

    def computeRdf(hdf, rdfGrid = None, saveRdf = True):
        """ Computes the radial distribution function and returns it as a numpy array.

        Positional arguments:  
        hdf -- Name of the hdf file
        
        Keyword arguments:
        rdfGrid -- A numpy array specifing the grid used to compute the radial distribution function
        (default None, i.e. rdfGrid = np.arange(0.01, 10.0, 0.01)  
        saveRdf -- Flag marking whether the rdf data should be added to the hdf file (default = True)
        
        """
        wvg = pd.read_hdf(hdf, "wvg")[0].to_numpy()
        ssf = pd.read_hdf(hdf, "ssf")[0].to_numpy()
        if (rdfGrid is None) : rdfGrid = np.arange(0.01, 10.0, 0.01)
        rdf = qp.computeRdf(rdfGrid, wvg, ssf)
        if (saveRdf):
            pd.DataFrame(rdfGrid).to_hdf(hdf, key="rdfGrid", mode="r+")
            pd.DataFrame(rdf).to_hdf(hdf, key="rdf", mode="r+")
        return rdf

    def computeInternalEnergy(hdf):
        """ Computes the internal energy and returns it to output.

        Positional arguments:  
        hdf -- Name of the hdf file
        
        """
        wvg = pd.read_hdf(hdf, "wvg")[0].to_numpy()
        ssf = pd.read_hdf(hdf, "ssf")[0].to_numpy()
        coupling = pd.read_hdf(hdf, "inputs")["coupling"][0].tolist()
        return qp.computeInternalEnergy(wvg, ssf, coupling)

    
class Plot():

    """ Class to gather methods used for plotting """
    
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
