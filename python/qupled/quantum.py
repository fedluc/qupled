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
import qupled.classic as classic

class Qstls(classic.Stls):

    """Class to solve the QSTLS scheme.

    Class used to setup and solve the quantum QSTLS scheme as described by
    `Schweng and Bohm <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037>`_ 
    This class inherits most of its methods and attributes from :obj:`~qupled.classic.Stls`

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
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1,
                 scheme2DIntegrals : str = "full",
                 threads : int = 1):
        # Allowed theories
        self.allowedTheories = ["QSTLS"]
        # Set theory
        self.inputs : qupled.qupled.QstlsInput = qp.QstlsInput() #: Inputs to solve the scheme.
        super()._setInputs(coupling, degeneracy, "QSTLS", chemicalPotential,
                           cutoff, error, mixing, iterations, matsubara,
                           outputFrequency, recoveryFile, resolution)
        self.inputs.int2DScheme = scheme2DIntegrals
        self.inputs.threads = threads
        if (fixed is not None): self.inputs.fixed = fixed
        # if (guess is not None): self.inputs.guess = guess
        # File to store output on disk
        self.hdfFileName = None

    # Compute
    def compute(self) -> None:
        """ Solves the scheme and saves the results to and hdf file. See the method :func:`~qupled.quantum.Qstls.save`
        to see which results are saved
        """
        self._checkInputs()
        self._unpackFixedAdrFiles()
        self.scheme = qp.Qstls(self.inputs)
        status = self.scheme.compute()
        self._checkStatusAndClean(status)
        self._setHdfFile()
        self._save()

    # Unpack zip folder with fixed component of the auxiliary density response
    # This is only a hook to the corresponding method in QstlsIet
    def _unpackFixedAdrFiles(self) -> None:
        pass
    
    # Save results to disk
    def _save(self) -> None:
        """ Stores the results obtained by solving the scheme. Extends :func:`~qupled.classic.Stls.save`
        by adding the option to save the auxiliary density response as a new dataframe in the hdf file. The
        auxiliary density response dataframe can be accessed as `adr`
        """
        super()._save()
        pd.DataFrame(self.scheme.adr).to_hdf(self.hdfFileName, key="adr")
        
    # Plot results
    def plot(self, toPlot : list[str],  matsubara : np.ndarray = None, rdfGrid : np.ndarray= None) -> None:
        """ Plots the results obtained stored in :obj:`~qupled.classic.Stls.scheme`. Extends
        :func:`~qupled.classic.Stls.plot` by adding the option to plot the auxiliary density
        response by passing `adr` to toPlot
        """
        super().plot(toPlot, matsubara, rdfGrid)
        if ("adr" in toPlot): self._plotAdr(matsubara)

    def _plotAdr(self, matsubara : list[int]) -> None:
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
    its methods and attributes from :obj:`~qupled.quantum.Qstls`

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
                         fixed, mixing, iterations,
                         matsubara, outputFrequency,
                         recoveryFile, resolution)
        # Allowed theories
        self.allowedTheories = ["QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"]
        # Set theory
        self.inputs : qupled.qupled.QstlsInput = qp.QstlsInput() #: Inputs to solve the scheme.
        super()._setInputs(coupling, degeneracy, theory, chemicalPotential,
                           cutoff, error, mixing, iterations, matsubara,
                           outputFrequency, recoveryFile, resolution)
        self.inputs.int2DScheme = scheme2DIntegrals
        self.inputs.threads = threads
        if (fixediet is not None): self.inputs.fixediet = fixediet
        self.inputs.iet = mapping
        self._checkInputs()
        # Temporary folder to store the unpacked files with the auxiliary density response
        self.fixediet = None
        self.tmpRunDir = None
        # File to store output on disk
        self.hdfFileName = None


    # Unpack zip folder with fixed component of the auxiliary density response
    def _unpackFixedAdrFiles(self) -> None:
        """ Unpacks the zip file storing the fixed component of the auxiliary density response """
        if (self.qInputs.fixediet != ""):
            self.tmpRunDir = "qupled_tmp_run_directory"
            zipFile = zf.ZipFile(self.qInputs.fixediet, "r")
            zipFile.extractall(self.tmpRunDir)
            self.qInputs.fixediet = self.tmpRunDir
    
    # Check that the dielectric scheme was solved without errors
    def _checkStatusAndClean(self, status) -> None:
        if (self.fixediet is not None):
            rmtree(self.tmpRunDir)
        if (status == 0):
            if os.path.isfile(self.scheme.recovery) : os.remove(self.scheme.recovery)
            print("Dielectric theory solved successfully!")
        else:
            sys.exit("Error while solving the dielectric theory")

            
    # Save results to disk
    def _save(self) -> None:
        """ Stores the results obtained by solving the scheme. Extends the corresponding method in the parent class
        by:  
        adding the option to save the bridge function adder as a new dataframe in the hdf file which can be
        accessed as bf  
        creating a zip file to group all the files produced at run-time and containing the fixed component of
        the auxiliary density response

        Stores the results obtained by solving the scheme. Extends :func:`~qupled.quantum.Qstls.save`
        by adding two functionalities: (1) save the bridge function adder as a new dataframe in the hdf file. The
        bridge function adder dataframe can be accessed as `bf` (2) create a zip file to group all the files
        produced at run-time and containing the fixed component of the auxiliary density response for the
        IET schemes
        
        """
        super()._save()
        pd.DataFrame(self.scheme.bf).to_hdf(self.hdfFileName, key="bf")
        # Zip all files for the fixed component of the auxiliary density response
        if (self.schemeqInputs.fixediet == ""):
            adrFileName = "adr_fixed_rs%5.3f_theta%5.3f_%s" % (self.inputs.coupling,
                                                               self.inputs.degeneracy,
                                                               self.inputs.theory)
            zipFile = zf.ZipFile(adrFileName + ".zip", "w")
            for adrFile in glob(adrFileName + "_wv*.bin"):
                zipFile.write(adrFile)
                os.remove(adrFile)

                
    # Plot results        
    def plot(self, toPlot, matsubara : list[int] = None, rdfGrid : np.ndarray= None) -> None:
        """ Plots the results obtained stored in :obj:`~qupled.classic.Stls.scheme`. Extends 
        :func:`~qupled.quantum.Qstls.plot` by adding the option to plot the bridge function
        adder by passing `bf` to toPlot
        """
        super().plot(toPlot, matsubara, rdfGrid)
        if ("bf" in toPlot):
            Plot.plot1D(self.scheme.wvg, self.scheme.bf, "Wave vector", "Bridge function adder")
        
