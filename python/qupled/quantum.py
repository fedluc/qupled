import sys
import os
from shutil import rmtree
from glob import glob
import zipfile as zf
import numpy as np
import pandas as pd
import qupled.util as qu
import qupled.qupled as qp
import qupled.classic as qc

# -----------------------------------------------------------------------
# Qstls class
# -----------------------------------------------------------------------

class Qstls(qc.Stls):

    """

    Class used to setup and solve the quantum QSTLS scheme as described by
    `Schweng and Bohm <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037>`_ 
    This class inherits most of its methods and attributes from :obj:`qupled.classic.Stls`.

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.  
        chemicalPotential: Initial guess for the chemical potential, defaults to [-10.0, 10.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for covergence, defaults to 1.0e-5.
        fixed: The name of the file storing the fixed component of the auxiliary density response.
               if no name is given the fixed component if computed from scratch
        mixing: Mixing parameter for iterative solution, defaults to 1.0.  
        guess:  Initial guess for the iterative solution, defaults to None, i.e. ssf from stls solution.
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
        scheme2DIntegrals: numerical scheme used to solve two-dimensional integrals. See :func:`qupled.qupled.Input.int2DScheme`
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
                 guess : qp.QstlsGuess = None,
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
        self._setInputs(coupling, degeneracy, "QSTLS", chemicalPotential, cutoff,
                        error, fixed, mixing, guess, iterations, matsubara,
                        outputFrequency, recoveryFile, resolution, scheme2DIntegrals,
                        threads)
        # Scheme to solve
        self.scheme : qp.Qstls = None
        # File to store output on disk
        self.hdfFileName = None

    # Setup inputs object
    def _setInputs(self,
                   coupling : float,
                   degeneracy : float,
                   theory : str,
                   chemicalPotential : list[float],
                   cutoff : float,
                   error : float,
                   fixed : str,
                   mixing : float,
                   guess : qp.QstlsGuess,
                   iterations : int,
                   matsubara : int,
                   outputFrequency : int,
                   recoveryFile : str,
                   resolution : float,
                   scheme2DIntegrals : str,
                   threads : int) -> None:
        super()._setInputs(coupling, degeneracy, theory, chemicalPotential,
                           cutoff, error, mixing, None, iterations, matsubara,
                           outputFrequency, recoveryFile, resolution)
        self.inputs.int2DScheme = scheme2DIntegrals
        self.inputs.threads = threads
        if (guess is not None): self.inputs.guess = guess
        if (fixed is not None): self.inputs.fixed = fixed
        
    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """ Solves the scheme and saves the results to an hdf file. Extends the output produced by 
        :func:`qupled.classic.Stls.compute` by adding the option to save the auxiliary density response
        as a new dataframe in the hdf file. The auxiliary density response dataframe can be accessed as `adr`.
        """
        self._checkInputs()
        self.scheme = qp.Qstls(self.inputs)
        status = self.scheme.compute()
        self._checkStatusAndClean(status)
        self._setHdfFile()
        self._save()

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """ Stores the results obtained by solving the scheme. 
        """
        super()._save()
        pd.DataFrame(self.scheme.adr).to_hdf(self.hdfFileName, key="adr")

    # Set the initial guess from a dataframe produced in output
    def setGuess(self, fileName : str) -> None:
        guess = qp.QstlsGuess()
        hdfData = qu.Hdf().read(fileName, ["wvg", "ssf"])
        guess.wvg = hdfData["wvg"]
        guess.ssf = hdfData["ssf"]
        self.inputs.guess = guess
        

# -----------------------------------------------------------------------
# QstlsIet class
# -----------------------------------------------------------------------
        
class QstlsIet(Qstls):

    """
    Class used to setup and solve the classical STLS-IET scheme as described by
    `Tolias <https://pubs.aip.org/aip/jcp/article/158/14/141102/
    2877795/Quantum-version-of-the-integral-equation-theory>`_. This class inherits most of
    its methods and attributes from :obj:`qupled.quantum.Qstls`.

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
        mapping: Classical to quantum mapping. See :func:`qupled.qupled.StlsInput.iet`
        mixing: Mixing parameter for iterative solution, defaults to 1.0.  
        guess:  Initial guess for the iterative solution, defaults to None, i.e. ssf from stls solution.
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
        scheme2DIntegrals: numerical scheme used to solve two-dimensional integrals. See :func:`qupled.qupled.Input.int2DScheme`
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
                 guess : qp.QstlsGuess = None,
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1,
                 scheme2DIntegrals : str = "full",
                 threads : int = 1):
        # Allowed theories
        self.allowedTheories = ["QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"]
        # Set theory
        self.inputs : qupled.qupled.QstlsInput = qp.QstlsInput() #: Inputs to solve the scheme.
        self._setInputs(coupling, degeneracy, theory, chemicalPotential, cutoff,
                        error, fixed, mixing, guess, iterations, matsubara,
                        outputFrequency, recoveryFile, resolution, scheme2DIntegrals,
                        threads)
        if (fixediet is not None): self.inputs.fixediet = fixediet
        self.inputs.iet = mapping
        self._checkInputs()
        # Temporary folder to store the unpacked files with the auxiliary density response
        self.fixediet = None
        self.tmpRunDir = None
        # Scheme to solve
        self.scheme : qp.QstlsIet = None
        # File to store output on disk
        self.hdfFileName = None


    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """ Solves the scheme and saves the results to an hdf file. Extends the output produced by 
        :func:`qupled.classic.Qstls.compute` by adding  by adding two functionalities: (1) save the
        bridge function adder as a new dataframe in the hdf file. The bridge function adder dataframe
        can be accessed as `bf` (2) create a zip file to group all the files produced at run-time
        and containing the fixed component of the auxiliary density response for the IET schemes.
        """
        self._checkInputs()
        self._setFixedIetFileName()
        self.scheme = qp.Qstls(self.inputs)
        status = self.scheme.compute()
        self._checkStatusAndClean(status)
        self._setHdfFile()
        self._save()

    # Set name of the file with the fixed component of the auxiliary density response
    @qu.MPI.synchronizeRanks
    def _setFixedIetFileName(self) -> None:
        """ Sets the file name for the file storing the fixed component of the auxiliary density response """
        if (self.inputs.fixediet != ""):
            self.tmpRunDir = "qupled_tmp_run_directory"
            self._unpackFixedAdrFiles()
            self.inputs.fixediet = self.tmpRunDir
        
    # Unpack zip folder with fixed component of the auxiliary density response
    @qu.MPI.runOnlyOnRoot
    def _unpackFixedAdrFiles(self) -> None:
        """ Unpacks the zip file storing the fixed component of the auxiliary density response """
        assert(self.inputs.fixediet != "")
        assert(self.tmpRunDir is not None)
        with zf.ZipFile(self.inputs.fixediet, "r") as zipFile:
            zipFile.extractall(self.tmpRunDir)
        
    # Check that the dielectric scheme was solved without errors
    @qu.MPI.runOnlyOnRoot
    def _checkStatusAndClean(self, status) -> None:
        # Remove the temporary run directory
        if (self.tmpRunDir is not None and os.path.isdir(self.tmpRunDir)) :
            rmtree(self.tmpRunDir)
        # Check that the scheme was solved correctly
        if (status == 0):
            if os.path.isfile(self.scheme.recovery) : os.remove(self.scheme.recovery)
            print("Dielectric theory solved successfully!")
        else:
            sys.exit("Error while solving the dielectric theory")

            
    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """ 
        Stores the results obtained by solving the scheme. 
        """
        super()._save()
        pd.DataFrame(self.scheme.bf).to_hdf(self.hdfFileName, key="bf")
        # Zip all files for the fixed component of the auxiliary density response
        if (self.inputs.fixediet == ""):
            adrFileName = "adr_fixed_rs%5.3f_theta%5.3f_%s" % (self.inputs.coupling,
                                                               self.inputs.degeneracy,
                                                               self.inputs.theory)
            with zf.ZipFile(adrFileName + ".zip", "w") as zipFile:
                for adrFile in glob(adrFileName + "_wv*.bin"):
                    zipFile.write(adrFile)
                    os.remove(adrFile)

    # Set the initial guess from a dataframe produced in output
    def setGuess(self, fileName : str) -> None:
        guess = qp.QstlsGuess()
        hdfData = qu.Hdf().read(fileName, ["wvg", "ssf", "adr", "matsubara"])
        guess.wvg = hdfData["wvg"]
        guess.ssf = hdfData["ssf"]
        guess.adr = np.ascontiguousarray(hdfData["adr"])
        guess.matsubara = hdfData["matsubara"]
        self.inputs.guess = guess
        
class QVSStls(qc.VSStls, Qstls):

    """
    Class used to setup and solve the quantum VS-STLS scheme as described by
    `Hayashi and Shimizu for ground state calculations <https://journals.jps.jp/doi/pdf/10.1143/JPSJ.48.16>`_ 
    and by ---Placeholder for the finite temperature qVS paper---.
    This class inherits most of its methods and attributes from :obj:`qupled.quantum.Qstls`.

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.  
        chemicalPotential: Initial guess for the chemical potential, defaults to [-10.0, 10.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for covergence, defaults to 1.0e-5.
        fixed: The name of the file storing the fixed component of the auxiliary density response.
               if no name is given the fixed component if computed from scratch
        mixing: Mixing parameter for iterative solution, defaults to 1.0.  
        guess:  Initial guess for the iterative solution, defaults to None, i.e. ssf from stls solution.
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
        scheme2DIntegrals: numerical scheme used to solve two-dimensional integrals. See :func:`qupled.qupled.Input.int2DScheme`
        alpha: Initial guess for the free parameter, defaults to [0.5, 1.0]
        couplingResolution: Resolution of the coupling parameter grid, defaults to 0.01
        degeneracyResolution: Resolution of the degeneracy parameter grid, defaults to 0.01
        errorAlpha: Minimum error for convergence in the free parameter iterations, defaults to 1.0e-3
        iterationsAlpha: Maximum number of iterations for the free parameter, defaults to 50
        errorIntegrals: Accuracy (as a relative error) for the integral computations, defaults to 1.0-5
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
                 guess : qp.QstlsGuess = None,
                 iterations : int = 1000,
                 matsubara : int = 128,
                 outputFrequency : int = 10,
                 recoveryFile : str = None,
                 resolution : float = 0.1,
                 scheme2DIntegrals : str = "full",
                 alpha : list[float] = [0.5,1.0],
                 couplingResolution : float = 0.1,
                 degeneracyResolution : float = 0.1,
                 errorAlpha : float = 1.0e-3,
                 iterationsAlpha : int = 50,
                 errorIntegrals : float = 1.0e-5,
                 threads : int = 1):
        # Allowed theories
        self.allowedTheories = ["QVSSTLS"]
        # Set theory
        self.inputs : qupled.qupled.QVSStlsInput = qp.QVSStlsInput() #: Inputs to solve the scheme.
        self._setInputs(coupling, degeneracy, "QVSSTLS", chemicalPotential, cutoff,
                        error, fixed, mixing, guess, iterations, matsubara,
                        outputFrequency, recoveryFile, resolution, scheme2DIntegrals,
                        alpha, couplingResolution, degeneracyResolution, errorAlpha,
                        iterationsAlpha, errorIntegrals, threads)
        # Scheme to solve
        self.scheme : qp.QVSStls = None
        # File to store output on disk
        self.hdfFileName = None

        # Setup inputs object
    def _setInputs(self,
                   coupling : float,
                   degeneracy : float,
                   theory : str,
                   chemicalPotential : list[float],
                   cutoff : float,
                   error : float,
                   fixed : str,
                   mixing : float,
                   guess : qp.QstlsGuess,
                   iterations : int,
                   matsubara : int,
                   outputFrequency : int,
                   recoveryFile : str,
                   resolution : float,
                   scheme2DIntegrals : str,
                   alpha : list[float],
                   couplingResolution : float,
                   degeneracyResolution : float,
                   errorAlpha : float,
                   iterationsAlpha : int,
                   errorIntegrals : float,
                   threads : int) -> None:
        super()._setInputs(coupling, degeneracy, theory, chemicalPotential,
                   cutoff, error, fixed, mixing, guess, iterations, matsubara,
                   outputFrequency, recoveryFile, resolution, scheme2DIntegrals, threads)
        self.inputs.alpha = alpha
        self.inputs.couplingResolution = couplingResolution
        self.inputs.degeneracyResolution = degeneracyResolution
        self.inputs.errorAlpha = errorAlpha
        self.inputs.iterationsAlpha = iterationsAlpha
        self.inputs.intError = errorIntegrals

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """ Solves the scheme and saves the results to an hdf file. Exactly as it is done by 
        :func:`qupled.classic.QStls.compute` .The auxiliary density response dataframe can be accessed as `adr`.
        """
        self._checkInputs()
        self.scheme = qp.QVSStls(self.inputs)
        status = self.scheme.compute()
        self._checkStatusAndClean(status)
        self._setHdfFile()
        self._save()

   