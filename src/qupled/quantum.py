from __future__ import annotations

import sys
import os
import shutil
import glob
from zipfile import ZipFile
import numpy as np
import pandas as pd
from qupled import native
import qupled.util as qu
import qupled.classic as qc

# -----------------------------------------------------------------------
# _QuantumIterativeScheme class
# -----------------------------------------------------------------------


class _QuantumIterativeScheme(qc._IterativeScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def getInitialGuess(fileName: str) -> _QuantumIterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            fileName : name of the file used to extract the information for the initial guess.
        """
        hdfData = qu.HDF().read(
            fileName,
            [
                qu.HDF.EntryKeys.WVG.value,
                qu.HDF.EntryKeys.SSF.value,
                qu.HDF.EntryKeys.ADR.value,
                qu.HDF.EntryKeys.MATSUBARA.value,
            ],
        )
        return _QuantumIterativeScheme.Guess(
            hdfData[qu.HDF.EntryKeys.WVG.value],
            hdfData[qu.HDF.EntryKeys.SSF.value],
            np.ascontiguousarray(hdfData[qu.HDF.EntryKeys.ADR.value]),
            hdfData[qu.HDF.EntryKeys.MATSUBARA.value],
        )

    # Save results to disk
    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        if scheme.inputs.degeneracy > 0:
            pd.DataFrame(scheme.adr).to_hdf(
                self.hdfFileName, key=qu.HDF.EntryKeys.ADR.value
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

        def toNative(self) -> native.QStlsGuess:
            native_guess = native.QstlsGuess()
            for attr, value in self.__dict__.items():
                native_value = value if value is not None else np.empty(0)
                setattr(native_guess, attr, native_value)
            return native_guess


# -----------------------------------------------------------------------
# Qstls class
# -----------------------------------------------------------------------


class Qstls(_QuantumIterativeScheme):

    # Compute
    @qu.MPI.record_time
    @qu.MPI.synchronize_ranks
    def compute(self, inputs: Qstls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.Qstls(inputs.toNative())
        self._compute(scheme)
        self._save(scheme)

    # Input class
    class Input(qc.Stls.Input):
        """
        Class used to manage the input for the :obj:`qupled.quantum.Qstls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.fixed: str = ""
            """ Name of the file storing the fixed component of the auxiliary density 
	    response in the QSTLS scheme. """
            self.guess: Qstls.Guess = Qstls.Guess()
            """Initial guess. Default = ``Qstls.Guess()``"""
            # Undocumented default values
            self.theory = "QSTLS"

        def toNative(self) -> native.QstlsInput:
            native_input = native.QstlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.toNative())
                else:
                    setattr(native_input, attr, value)
            return native_input


# -----------------------------------------------------------------------
# QstlsIet class
# -----------------------------------------------------------------------


class QstlsIet(_QuantumIterativeScheme):
    """
    Args:
        inputs: Input parameters.
    """

    # Compute
    @qu.MPI.record_time
    @qu.MPI.synchronize_ranks
    def compute(self, inputs: QstlsIet.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self._unpackFixedAdrFiles(inputs)
        scheme = native.Qstls(inputs.toNative())
        self._compute(scheme)
        self._save(scheme)
        self._zipFixedAdrFiles(inputs)
        self._cleanFixedAdrFiles(scheme.inputs)

    # Unpack zip folder with fixed component of the auxiliary density response
    @qu.MPI.run_only_on_root
    def _unpackFixedAdrFiles(self, inputs) -> None:
        fixedIetSourceFile = inputs.fixed_iet
        if inputs.fixed_iet != "":
            inputs.fixed_iet = "qupled_tmp_run_directory"
        if fixedIetSourceFile != "":
            with ZipFile(fixedIetSourceFile, "r") as zipFile:
                zipFile.extractall(inputs.fixed_iet)

    # Save results to disk
    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        super()._save(scheme)
        pd.DataFrame(scheme.bf).to_hdf(self.hdfFileName, key="bf")

    # Zip all files for the fixed component of the auxiliary density response
    @qu.MPI.run_only_on_root
    def _zipFixedAdrFiles(self, inputs) -> None:
        if inputs.fixed_iet == "":
            degeneracy = inputs.degeneracy
            matsubara = inputs.matsubara
            theory = inputs.theory
            adrFile = f"adr_fixed_theta{degeneracy:5.3f}_matsubara{matsubara}_{theory}"
            adrFileZip = f"{adrFile}.zip"
            adrFileBin = f"{adrFile}_wv*.bin"
            with ZipFile(adrFileZip, "w") as zipFile:
                for binFile in glob.glob(adrFileBin):
                    zipFile.write(binFile)
                    os.remove(binFile)

    # Remove temporaray run directory
    @qu.MPI.run_only_on_root
    def _cleanFixedAdrFiles(self, inputs) -> None:
        if os.path.isdir(inputs.fixed_iet):
            shutil.rmtree(inputs.fixed_iet)

    # Input class
    class Input(qc.StlsIet.Input, Qstls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.QStlsIet` class.
        Accepted theories: ``QSTLS-HNC``, ``QSTLS-IOI`` and ``QSTLS-LCT``.
        """

        def __init__(self, coupling: float, degeneracy: float, theory: str):
            qc.StlsIet.Input.__init__(self, coupling, degeneracy, "STLS-HNC")
            Qstls.Input.__init__(self, coupling, degeneracy)
            if theory not in {"QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"}:
                sys.exit("Invalid dielectric theory")
            self.theory = theory
            self.fixed_iet = ""
            """
            Name of the zip file storing the iet part of the fixed components
            of the auxiliary density response. Default = ``""``
            """

        def toNative(self) -> native.QstlsInput:
            native_input = native.QstlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.toNative())
                else:
                    setattr(native_input, attr, value)
            return native_input


# -----------------------------------------------------------------------
# QVSStls class
# -----------------------------------------------------------------------


class QVSStls(_QuantumIterativeScheme):

    # Compute
    @qu.MPI.record_time
    @qu.MPI.synchronize_ranks
    def compute(self, inputs: QVSStls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self._unpackFixedAdrFiles(inputs)
        scheme = native.QVSStls(inputs.toNative())
        self._compute(scheme)
        self._save(scheme)
        self._zipFixedAdrFiles(inputs)
        self._cleanFixedAdrFiles(scheme.inputs)

    # Unpack zip folder with fixed component of the auxiliary density response
    @qu.MPI.run_only_on_root
    def _unpackFixedAdrFiles(self, inputs) -> None:
        fixedSourceFile = inputs.fixed
        if inputs.fixed != "":
            inputs.fixed = "qupled_tmp_run_directory"
        if fixedSourceFile != "":
            with ZipFile(fixedSourceFile, "r") as zipFile:
                zipFile.extractall(inputs.fixed)

    # Save results to disk
    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        super()._save(scheme)
        pd.DataFrame(scheme.free_energy_grid).to_hdf(
            self.hdfFileName, key=qu.HDF.EntryKeys.FXC_GRID.value
        )
        pd.DataFrame(scheme.free_energy_integrand).to_hdf(
            self.hdfFileName, key=qu.HDF.EntryKeys.FXCI.value
        )
        pd.DataFrame(scheme.alpha).to_hdf(
            self.hdfFileName, key=qu.HDF.EntryKeys.ALPHA.value
        )

    # Zip all files for the fixed component of the auxiliary density response
    @qu.MPI.run_only_on_root
    def _zipFixedAdrFiles(self, inputs) -> None:
        if inputs.fixed == "":
            degeneracy = inputs.degeneracy
            matsubara = inputs.matsubara
            theory = inputs.theory
            adrFileZip = (
                f"adr_fixed_theta{degeneracy:5.3f}_matsubara{matsubara}_{theory}.zip"
            )
            adrFileBin = "THETA*.bin"
            with ZipFile(adrFileZip, "w") as zipFile:
                for binFile in glob.glob(adrFileBin):
                    zipFile.write(binFile)
                    os.remove(binFile)

    # Remove the temporary run directory
    @qu.MPI.run_only_on_root
    def _cleanFixedAdrFiles(self, inputs) -> None:
        if os.path.isdir(inputs.fixed):
            shutil.rmtree(inputs.fixed)

    # Set the free energy integrand from a dataframe produced in output
    @staticmethod
    def getFreeEnergyIntegrand(fileName: str) -> native.FreeEnergyIntegrand:
        return qc.VSStls.getFreeEnergyIntegrand(fileName)

    # Input class
    class Input(qc.VSStls.Input, Qstls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.QVSStls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            qc.VSStls.Input.__init__(self, coupling, degeneracy)
            Qstls.Input.__init__(self, coupling, degeneracy)
            # Undocumented default values
            self.theory: str = "QVSSTLS"

        def toNative(self) -> native.QVSStlsInput:
            native_input = native.QVSStlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.toNative())
                else:
                    setattr(native_input, attr, value)
            return native_input
