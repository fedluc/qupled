#!/usr/bin/env python

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
        error: Minimum error for covergence, defaults to 1.0e-5.
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
                 guess : qp.SlfcGuess = None,
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1 ):
        # Allowed theories
        self.allowedTheories : list[str] = ["STLS"]
        # Input object
        self.inputs : qp.StlsInput = qp.StlsInput(coupling, degeneracy, "STLS") #: Inputs to solve the scheme.
        # Scheme to solve and associated input and solution
        self.scheme : qp.Stls = None #: Object that represents the scheme, performs the calculations and stores the solution.
        self.schemeInputs : qp.StlsInput = None
        # File to store output on disk
        self.hdfFileName : str = None #: Name of the hdf output file.
        # Optional parameters
        self.inputs.chemicalPotential = chemicalPotential
        self.inputs.cutoff = cutoff
        self.inputs.error = error
        if (guess is not None): self.inputs.guess = guess
        self.inputs.mixing = mixing
        self.inputs.iterations = iterations
        self.inputs.matsubara = matsubara
        self.inputs.outputFrequency = outputFrequency
        if (recoveryFile is not None): self.inputs.recoveryFile = recoveryFile
        self.inputs.resolution = resolution

    # Check input before computing
    def checkInputs(self) -> None:
        """ Checks that the content of :obj:`inputs` is correct """
        if (self.inputs.theory not in self.allowedTheories):
            sys.exit("Invalid dielectric theory")

    # Compute
    def compute(self) -> None:
        """ Solves the scheme and saves the results to and hdf file. See the method :func:`~qupled.Static.Stls.save`
        to see which results are saved
        """
        self.checkInputs()
        self.schemeInputs = self.inputs
        self.scheme = qp.Stls(self.schemeInputs)
        status = self.scheme.compute()
        self.checkStatusAndClean(status)        
        self.setHdfFile()
        self.save()

    # Check that the dielectric scheme was solved without errors
    def checkStatusAndClean(self, status : bool) -> None:
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
    def setHdfFile(self) -> None:
        """ Sets the name of the hdf file used to store the output """
        self.hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (self.schemeInputs.coupling,
                                                         self.schemeInputs.degeneracy,
                                                         self.schemeInputs.theory)
    
    def save(self) -> None:
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
            "coupling" : self.schemeInputs.coupling,
            "degeneracy" : self.schemeInputs.degeneracy,
            "theory" : self.schemeInputs.theory,
            "error" : self.schemeInputs.error,
            "resolution" : self.schemeInputs.resolution,
            "cutoff" : self.schemeInputs.cutoff,
            "matsubara" : self.schemeInputs.matsubara
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
        if (self.schemeInputs == None):
            sys.exit("No solution to compute the radial distribution function")
        rdf = self.scheme.getRdf(rdfGrid)
        if (writeToHdf) :
            pd.DataFrame(rdfGrid).to_hdf(self.hdfFileName, key="rdfGrid", mode="r+")
            pd.DataFrame(rdf).to_hdf(self.hdfFileName, key="rdf", mode="r+")
        return rdf
        
    # Plot results
    def plot(self, toPlot : list[str], matsubara : np.ndarray = None, rdfGrid : np.ndarray = None) -> None:
        """ Plots the results obtained stored in :obj:`scheme`.

        Args:  
            toPlot: A list of quantities to plot. This list can include idr (ideal density response), rdf
                (radial distribution function), sdr (static density response), slfc (static local field correction)
                ssf (static structure factor) and ssfHF (Hartree-Fock static structure factor)  
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Default = None, see :func:`~qupled.Static.Stls.plotIdr`)  
            rdfGrid: The grid used to compute the radial distribution function. Applies only when the radial
                distribution function is plotted (Default = None, see :func:`~qupled.Static.Stls.computeRdf`)
        
        """
        wvg = self.scheme.wvg
        xlabel = "Wave vector"
        if ("idr" in toPlot):
            self.plotIdr(matsubara)
        if ("rdf" in toPlot):
            self.plotRdf(rdfGrid)
        if ("sdr" in toPlot):
            Plot.plot1D(wvg, self.scheme.sdr, xlabel, "Static density response")
        if ("slfc" in toPlot):
            Plot.plot1D(wvg, self.scheme.slfc, xlabel, "Static local field correction")
        if ("ssf" in toPlot):
            Plot.plot1D(wvg, self.scheme.ssf, xlabel, "Static structure factor")
        if ("ssfHF" in toPlot):
            Plot.plot1D(wvg, self.scheme.ssfHF, xlabel, "Hartree-Fock static structure factor")
        
    def plotIdr(self, matsubara : np.ndarray = None) -> None:
        """ Plots the ideal density response.
        
        Args:  
            matsubara:  A list of matsubara frequencies to plot. (Default =  all matsubara frequencies are plotted)
        
        """
        if (matsubara is None) : matsubara = np.arange(self.inputs.matsubara)
        Plot.plot1DParametric(self.scheme.wvg, self.scheme.idr,
                              "Wave vector", "Ideal density response",
                              matsubara)
        
    def plotRdf(self, rdfGrid : np.ndarray = None) -> None:
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
    This class inherits most of its methods and attributes from :obj:`~qupled.Static.Stls`

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.  
        chemicalPotential: Initial guess for the chemical potential, defaults to [-10.0, 10.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for covergence, defaults to 1.0e-5.
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
                 guess : qp.SlfcGuess = None,
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1,
                 scheme2DIntegrals : str = "full"):
        # Call parent constructor
        super().__init__(coupling, degeneracy,
                         chemicalPotential, cutoff, error,
                         mixing, guess, iterations,
                         matsubara, outputFrequency,
                         recoveryFile, resolution)
        # Allowed theories
        self.allowedTheories = ["STLS-HNC", "STLS-IOI", "STLS-LCT"]
        # Set theory
        self.inputs.theory = theory
        self.checkInputs()
        # File to store output on disk
        self.hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (self.inputs.coupling,
                                                         self.inputs.degeneracy,
                                                         self.inputs.theory)
        # Non-default inputs
        self.inputs.iet = mapping
        self.inputs.int2DScheme = scheme2DIntegrals
            
    # Plot results
    def plot(self, toPlot, matsubara : list[int] = None, rdfGrid : np.ndarray= None) -> None:
        """ Plots the results obtained stored in :obj:`~qupled.Static.Stls.scheme`. Extends 
        :func:`~qupled.Static.Stls.plot` by adding the option to plot the bridge function
        adder by passing `bf` to toPlot
        """
        super().plot(toPlot, matsubara, rdfGrid)
        if ("bf" in toPlot):
            Plot.plot1D(self.scheme.wvg, self.scheme.bf, "Wave vector", "Bridge function adder")
        
    # Save results to disk
    def save(self) -> None:
        """ Stores the results obtained by solving the scheme. Extends :func:`~qupled.Static.Stls.save`
        by adding the option to save the bridge function adder as a new dataframe in the hdf file. The
        bridge function adder dataframe can be accessed as `bf`
        """
        super().save()
        pd.DataFrame(self.scheme.bf).to_hdf(self.hdfFileName, key="bf")



class Qstls(Stls):

    """Class to solve the STLS-IET schemes.

    Class used to setup and solve the quantum QSTLS scheme as described by
    `Schweng and Bohm <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037>`_ 
    This class inherits most of its methods and attributes from :obj:`~qupled.Static.Stls`

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.  
        chemicalPotential: Initial guess for the chemical potential, defaults to [-10.0, 10.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for covergence, defaults to 1.0e-5.
        fixed: The name of the file storing the fixed component of the auxiliary density response.
               if no name is given the fixed component if computed from scratch
        mixing: Mixing parameter for iterative solution, defaults to 1.0.  
        guess:  Initial guess for the iterative solution, defaults to None, i.e. slfc = 0.
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
        scheme2DIntegrals: numerical scheme used to solve two-dimensional integrals. See :func:`~qupled.qupled.Input.int2DScheme`
        threads: OMP threads for parallel calculations
    """
    
    # Constructor
    def __init__(self,
                 coupling : float, 
                 degeneracy : float,
                 chemicalPotential : list[float] = [-10.0,10.0],
                 cutoff : float = 10.0,
                 error : float = 1.0e-5,
                 fixed : str = None,
                 mixing : float = 1.0,
                 guess : qp.SlfcGuess = None,
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1,
                 scheme2DIntegrals : str = "full",
                 threads : int = 1):
        
        # Call parent constructor
        super().__init__(coupling, degeneracy,
                         chemicalPotential, cutoff, error,
                         mixing, None, iterations,
                         matsubara, outputFrequency,
                         recoveryFile, resolution)
        # Allowed theories
        self.allowedTheories = ["QSTLS"]
        # Set theory
        self.inputs.theory = "QSTLS"
        # Qstls inputs
        self.qInputs : qp.QstlsInput = qp.QstlsInput() #: Inputs to solve the quantum schemes.
        self.schemeqInputs = None;
        self.checkInputs()
        # File to store output on disk
        self.hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (self.inputs.coupling,
                                                         self.inputs.degeneracy,
                                                         self.inputs.theory)
        # Non-default inputs
        if (fixed is not None): self.qInputs.fixed = fixed
        if (guess is not None): self.qInputs.guess = guess
        self.inputs.int2DScheme = scheme2DIntegrals
        self.inputs.threads = threads

    # Compute
    def compute(self) -> None:
        self.checkInputs()
        self.unpackFixedAdrFiles()
        self.schemeInputs = self.inputs
        self.schemeqInputs = self.qInputs
        self.scheme = qp.Qstls(self.schemeInputs, self.schemeqInputs)
        status = self.scheme.compute()
        self.checkStatusAndClean(status)
        self.setHdfFile()
        self.save()

    # Unpack zip folder with fixed component of the auxiliary density response
    # This is only a hook to the corresponding method in QstlsIet
    def unpackFixedAdrFiles(self) -> None:
        pass
    
    # Save results to disk
    def save(self) -> None:
        """ Stores the results obtained by solving the scheme. Extends :func:`~qupled.Static.Stls.save`
        by adding the option to save the auxiliary density response as a new dataframe in the hdf file. The
        auxiliary density response dataframe can be accessed as `adr`
        """
        super().save()
        pd.DataFrame(self.scheme.adr).to_hdf(self.hdfFileName, key="adr")
        
    # Plot results
    def plot(self, toPlot : list[str],  matsubara : np.ndarray = None, rdfGrid : np.ndarray= None) -> None:
        """ Plots the results obtained stored in :obj:`~qupled.Static.Stls.scheme`. Extends
        :func:`~qupled.Static.Stls.plot` by adding the option to plot the auxiliary density
        response by passing `adr` to toPlot
        """
        super().plot(toPlot, matsubara, rdfGrid)
        if ("adr" in toPlot): self.plotAdr(matsubara)

    def plotAdr(self, matsubara : list[int]) -> None:
        """ Plots the auxiliary density response.
        
        Args:  
            matsubara:  A list of matsubara frequencies to plot. (Default =  all matsubara frequencies are plotted)
        
        """
        if (matsubara is None) : matsubara = np.arange(self.inputs.matsubara)
        Plot.plot1DParametric(self.scheme.wvg, self.scheme.adr,
                              "Wave vector", "Auxiliary density response",
                              matsubara)
        

class QstlsIet(Qstls):

    """Class to solve the QSTLS-IET schemes.

    Class used to setup and solve the classical STLS-IET scheme as described by
    `Tolias <https://pubs.aip.org/aip/jcp/article/158/14/141102/
    2877795/Quantum-version-of-the-integral-equation-theory>`_. This class inherits most of
    its methods and attributes from :obj:`~qupled.Static.Qstls`

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.  
        chemicalPotential: Initial guess for the chemical potential, defaults to [-10.0, 10.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for covergence, defaults to 1.0e-5.
        fixed: The name of the file storing the fixed component of the auxiliary density response.
               if no name is given the fixed component is computed from scratch.
        fixediet: The name of the zip file storing the files with the fixed component of the auxiliary
                  density response for the IET schemes. If no name is given the fixed component
                  is computed from scratch.
        mapping: Classical to quantum mapping. See :func:`~qupled.qupled.StlsInput.iet`
        mixing: Mixing parameter for iterative solution, defaults to 1.0.  
        guess:  Initial guess for the iterative solution, defaults to None, i.e. slfc = 0.
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
        scheme2DIntegrals: numerical scheme used to solve two-dimensional integrals. See :func:`~qupled.qupled.Input.int2DScheme`
        threads: OMP threads for parallel calculations
    """
    # Constructor
    def __init__(self,
                 coupling : float, 
                 degeneracy : float,
                 theory : str,
                 chemicalPotential : list[float] = [-10.0,10.0],
                 cutoff : float = 10.0,
                 error : float = 1.0e-5,
                 fixed : str = None,
                 fixediet : str = None,
                 mapping : str = "standard",
                 mixing : float = 1.0,
                 guess : qp.SlfcGuess = None,
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1,
                 scheme2DIntegrals : str = "full",
                 threads : int = 1):
        
        # Call parent constructor
        super().__init__(coupling, degeneracy,
                         chemicalPotential, cutoff, error,
                         fixed, mixing, guess, iterations,
                         matsubara, outputFrequency,
                         recoveryFile, resolution)
        # Allowed theories
        self.allowedTheories = ["QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"]
        # Set theory
        self.inputs.theory = theory
        self.checkInputs()
        # Temporary folder to store the unpacked files with the auxiliary density response
        self.fixediet = None
        self.tmpRunDir = None
        # File to store output on disk
        self.hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (self.inputs.coupling,
                                                           self.inputs.degeneracy,
                                                           self.inputs.theory)
        # Non-default inputs
        if (fixediet is not None): self.fixediet = fixediet
        self.inputs.iet = mapping

    # Unpack zip folder with fixed component of the auxiliary density response
    def unpackFixedAdrFiles(self) -> None:
        """ Unpacks the zip file storing the fixed component of the auxiliary density response """
        if (self.fixediet is not None):
            self.tmpRunDir = "qupled_tmp_run_directory"
            zipFile = zf.ZipFile(self.fixediet, "r")
            zipFile.extractall(self.tmpRunDir)
            self.qInputs.fixediet = self.tmpRunDir
    
    # Check that the dielectric scheme was solved without errors
    def checkStatusAndClean(self, status) -> None:
        if (self.fixediet is not None):
            rmtree(self.tmpRunDir)
        if (status == 0):
            if os.path.isfile(self.scheme.recovery) : os.remove(self.scheme.recovery)
            print("Dielectric theory solved successfully!")
        else:
            sys.exit("Error while solving the dielectric theory")

            
    # Save results to disk
    def save(self) -> None:
        """ Stores the results obtained by solving the scheme. Extends the corresponding method in the parent class
        by:  
        adding the option to save the bridge function adder as a new dataframe in the hdf file which can be
        accessed as bf  
        creating a zip file to group all the files produced at run-time and containing the fixed component of
        the auxiliary density response

        Stores the results obtained by solving the scheme. Extends :func:`~qupled.Static.Qstls.save`
        by adding two functionalities: (1) save the bridge function adder as a new dataframe in the hdf file. The
        bridge function adder dataframe can be accessed as `bf` (2) create a zip file to group all the files
        produced at run-time and containing the fixed component of the auxiliary density response for the
        IET schemes
        
        """
        super().save()
        pd.DataFrame(self.scheme.bf).to_hdf(self.hdfFileName, key="bf")
        # Zip all files for the fixed component of the auxiliary density response
        if (self.schemeqInputs.fixediet == ""):
            adrFileName = "adr_fixed_rs%5.3f_theta%5.3f_%s" % (self.schemeInputs.coupling,
                                                               self.schemeInputs.degeneracy,
                                                               self.schemeInputs.theory)
            zipFile = zf.ZipFile(adrFileName + ".zip", "w")
            for adrFile in glob(adrFileName + "_wv*.bin"):
                zipFile.write(adrFile)
                os.remove(adrFile)

                
    # Plot results        
    def plot(self, toPlot, matsubara : list[int] = None, rdfGrid : np.ndarray= None) -> None:
        """ Plots the results obtained stored in :obj:`~qupled.Static.Stls.scheme`. Extends 
        :func:`~qupled.Static.Qstls.plot` by adding the option to plot the bridge function
        adder by passing `bf` to toPlot
        """
        super().plot(toPlot, matsubara, rdfGrid)
        if ("bf" in toPlot):
            Plot.plot1D(self.scheme.wvg, self.scheme.bf, "Wave vector", "Bridge function adder")
        


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
