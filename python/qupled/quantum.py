from __future__ import annotations

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
# QuantumIterativeScheme class
# -----------------------------------------------------------------------


class QuantumIterativeScheme(qc.IterativeScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def getInitialGuess(fileName: str) -> qp.QstlsGuess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            fileName : name of the file used to extract the information for the initial guess.
        """
        guess = qp.QstlsGuess()
        hdfData = qu.Hdf().read(fileName, ["wvg", "ssf"])
        guess.wvg = hdfData["wvg"]
        guess.ssf = hdfData["ssf"]
        return guess

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save()
        pd.DataFrame(self.adr).to_hdf(self.hdfFileName, key="adr")


# -----------------------------------------------------------------------
# Qstls class
# -----------------------------------------------------------------------


class QstlsMetaclass(type(QuantumIterativeScheme), type(qp.Qstls)):
    pass


class Qstls(QuantumIterativeScheme, qp.Qstls, metaclass=QstlsMetaclass):
    """

    Class used to setup and solve the quantum QSTLS scheme as described by
    `Schweng and Bohm <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037>`_
    This class inherits most of its methods and attributes from :obj:`qupled.classic.Stls`.

    Args:
        inputs: Input parameters used to solve the scheme.
    """

    class Input(qc.Stls.Input, qp.QstlsInput):
        """
        Class used to manage the input for the :obj:`qupled.quantum.Qstls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy, False)
            self.fixed: str = ""
            """ Name of the file storing the fixed component of the auxiliary density 
	    response in the QSTLS scheme. """
            self.guess: qp.QstlsGuess = qp.QstlsGuess()
            """ Initial guess """
            # Undocumented default values
            self.theory = "QSTLS"

    # Constructor
    def __init__(self, inputs: Qstls.Input):
        # Construct the base classes
        super().__init__(inputs)
        # File to store output on disk
        self.hdfFileName: str = self._getHdfFile()  #: Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results.

        The results are stored as pandas dataframes in an hdf file with the following keywords:

        - info: A dataframe containing information on the input parameters, it includes:

          - coupling: the coupling parameter,
          - degeneracy: the degeneracy parameter,
          - error: the residual error at the end of the solution
          - theory: the theory that is being solved,
          - resolution: the resolution in the wave-vector grid,
          - cutoff: the cutoff in the wave-vector grid,
          - matsubara: the number of matsubara frequencies

        - adr (*ndarray*, 2D): the auxiliary density response
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

        The name of the hdf file is stored in :obj:`hdfFileName`.
        """
        super().computeScheme(super().compute, self._save)


# -----------------------------------------------------------------------
# QstlsIet class
# -----------------------------------------------------------------------


class QstlsIet(QuantumIterativeScheme, qp.Qstls, metaclass=QstlsMetaclass):
    """
    Class used to setup and solve the classical QSTLS-IET scheme as described by
    `Tolias <https://pubs.aip.org/aip/jcp/article/158/14/141102/
    2877795/Quantum-version-of-the-integral-equation-theory>`_.


    Args:
        inputs: Input parameters used to solve the scheme.
    """

    class Input(Qstls.Input, qp.QstlsInput):
        """
        Class used to manage the input for the :obj:`qupled.quantum.QStlsIet` class.
        """

        def __init__(self, coupling: float, degeneracy: float, theory: str):
            super().__init__(coupling, degeneracy)
            if theory not in {"QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"}:
                print(theory)
                sys.exit("Invalid dielectric theory")
            self.theory = theory
            """ Dielectric theory to solve  """
            self.mapping = "standard"
            """ Classical-to-quantum mapping used in the iet schemes
            allowed options include:
        
              - standard: inversely proportional to the degeneracy parameter
        
	      - sqrt: inversely proportional to the square root of the sum
                      of the squares of one and the degeneracy parameter
            
              - linear: inversely proportional to one plus the degeneracy
                        parameter.
            
            Far from the ground state all mappings lead identical results, but at
            the ground state they can differ significantly (the standard
            mapping diverges)
            """
            self.fixediet = ""
            """ Name of the zip file storing the fixed components of the auxiliary density
	    response in the QSTLS-IET schemes """

    # Constructor
    def __init__(self, inputs: QstlsIet.Input):
        # Setup the folder structure to load the fixed component of the adr
        self.fixedIetSourceFile = inputs.fixediet
        if inputs.fixediet != "":
            inputs.fixediet = "qupled_tmp_run_directory"
        # Construct the base class
        super().__init__(inputs)
        # File to store the output on disk
        self.hdfFileName: str = self._getHdfFile()  # . Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results.

        The results are stored as pandas dataframes in an hdf file with the following keywords:

        - info: A dataframe containing information on the input parameters, it includes:

          - coupling: the coupling parameter,
          - degeneracy: the degeneracy parameter,
          - error: the residual error at the end of the solution
          - theory: the theory that is being solved,
          - resolution: the resolution in the wave-vector grid,
          - cutoff: the cutoff in the wave-vector grid,
          - matsubara: the number of matsubara frequencies

        - adr (*ndarray*, 2D): the auxiliary density response
        - bf (*ndarray*): the bridge function adder
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

        The name of the hdf file is stored in :obj:`hdfFileName`.
        """
        self._unpackFixedAdrFiles()
        status = super().compute()
        self._checkStatusAndClean(status)
        self._save()

    # Unpack zip folder with fixed component of the auxiliary density response
    @qu.MPI.runOnlyOnRoot
    def _unpackFixedAdrFiles(self) -> None:
        if self.fixedIetSourceFile != "":
            with zf.ZipFile(self.fixedIetSourceFile, "r") as zipFile:
                zipFile.extractall(self.inputs.fixediet)

    # Check that the dielectric scheme was solved without errors
    @qu.MPI.runOnlyOnRoot
    def _checkStatusAndClean(self, status) -> None:
        # Remove the temporary run directory
        if os.path.isdir(self.inputs.fixediet):
            rmtree(self.inputs.fixediet)
        # Check that the scheme was solved correctly
        super()._checkStatusAndClean(status)

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        super()._save()
        pd.DataFrame(self.bf).to_hdf(self.hdfFileName, key="bf")
        # Zip all files for the fixed component of the auxiliary density response
        if self.inputs.fixediet == "":
            adrFileName = "adr_fixed_theta%5.3f_matsubara%d_%s" % (
                self.inputs.degeneracy,
                self.inputs.matsubara,
                self.inputs.theory,
            )
            with zf.ZipFile(adrFileName + ".zip", "w") as zipFile:
                for adrFile in glob(adrFileName + "_wv*.bin"):
                    zipFile.write(adrFile)
                    os.remove(adrFile)

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def getInitialGuess(fileName: str) -> qp.QstlsInput:
        guess = QuantumIterativeScheme.getInitialGuess(fileName)
        hdfData = qu.Hdf().read(fileName, ["adr", "matsubara"])
        guess.adr = np.ascontiguousarray(hdfData["adr"])
        guess.matsubara = hdfData["matsubara"]
        return guess


# -----------------------------------------------------------------------
# QVSStls class
# -----------------------------------------------------------------------

class QVSStlsMetaclass(type(QuantumIterativeScheme), type(qp.QVSStls)):
    pass


class QVSStls(QuantumIterativeScheme, qp.QVSStls, metaclass=QVSStlsMetaclass):
    """
    Class used to setup and solve the quantum VS-STLS scheme.

    Args:
        inputs: Input parameters used to solve the scheme.
    """

    class Input(Qstls.Input, qp.QVSStlsInput):
        """
        Class used to manage the input for the :obj:`qupled.classic.Stls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            """ Name of the theory that is solved """
            self.alpha: list[float] = [0.5, 1.0]
            """ Initial guess for the free parameter """
            self.couplingResolution: float = 0.1
            """ Resolution of the coupling parameter grid """
            self.degeneracyResolution: float = 0.1
            """ Resolution of the degeneracy parameter grid """
            self.errorAlpha: float = 1.0e-3
            """ Minimum error for convergence in the free parameter """
            self.iterationsAlpha: int = 50
            """ Maximum number of iterations to determine the free parameter """
            self.freeEnergyIntegrand: qupled.FreeEnergyIntegrand = (
                qp.FreeEnergyIntegrand()
            )
            """ Pre-computed free energy integrand """
            # Undocumented default values
            self.threads = 9
            self.theory = "QVSSTLS"
            
    # Constructor
    def __init__(self, inputs: Qstls.Input):
        # Setup the folder structure to load the fixed component of the adr
        self.fixedSourceFile = inputs.fixed
        if inputs.fixed != "":
            inputs.fixed = "qupled_tmp_run_directory"
        # Construct the base classes
        super().__init__(inputs)
        # File to store output on disk
        self.hdfFileName: str = self._getHdfFile()  #: Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results.

        The results are stored as pandas dataframes in an hdf file with the following keywords:

        - info: A dataframe containing information on the input parameters, it includes:

          - coupling: the coupling parameter,
          - degeneracy: the degeneracy parameter,
          - error: the residual error at the end of the solution
          - theory: the theory that is being solved,
          - resolution: the resolution in the wave-vector grid,
          - cutoff: the cutoff in the wave-vector grid,
          - matsubara: the number of matsubara frequencies

        - adr (*ndarray*, 2D): the auxiliary density response
        - fxcGrid (*ndarray*): coupling parameter grid
        - fxci (*ndarray*): the free energy integrand
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

        The name of the hdf file is stored in :obj:`hdfFileName`.
        """
        self._unpackFixedAdrFiles()
        status = super().compute()
        self._checkStatusAndClean(status)
        self._save()

    # Unpack zip folder with fixed component of the auxiliary density response
    @qu.MPI.runOnlyOnRoot
    def _unpackFixedAdrFiles(self) -> None:
        if self.fixedSourceFile != "":
            with zf.ZipFile(self.inputs.fixed, "r") as zipFile:
                zipFile.extractall(self.inputs.fixediet)
            
    # Check that the dielectric scheme was solved without errors
    @qu.MPI.runOnlyOnRoot
    def _checkStatusAndClean(self, status) -> None:
        # Remove the temporary run directory
        if os.path.isdir(self.inputs.fixed):
            rmtree(self.inputs.fixed)
        # Check that the scheme was solved correctly
        super()._checkStatusAndClean(status)

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        super()._save()
        pd.DataFrame(self.freeEnergyGrid).to_hdf(self.hdfFileName, key="fxcGrid")
        pd.DataFrame(self.freeEnergyIntegrand).to_hdf(self.hdfFileName, key="fxci")
        pd.DataFrame(self.alpha).to_hdf(self.hdfFileName, key="alpha")
        # Zip all files for the fixed component of the auxiliary density response
        if self.inputs.fixed == "":
            adrFileName = "adr_fixed_theta%5.3f_matsubara%d.zip" % (
                self.inputs.degeneracy,
                self.inputs.matsubara,
            )
            with zf.ZipFile(adrFileName, "w") as zipFile:
                for adrFile in glob("THETA*.bin"):
                    zipFile.write(adrFile)
                    os.remove(adrFile)
                    
    # Set the free energy integrand from a dataframe produced in output
    @staticmethod
    def getFreeEnergyIntegrand(fileName: str) -> qp.FreeEnergyIntegrand():
        return qc.VSStls.getFreeEnergyIntegrand(fileName)
