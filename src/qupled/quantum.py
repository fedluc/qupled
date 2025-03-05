from __future__ import annotations

import glob
import os
import shutil
import sys
import zipfile

import numpy as np
import pandas as pd

import qupled.classic as qc
import qupled.util as qu
import qupled.native as qn

# -----------------------------------------------------------------------
# _QuantumIterativeScheme class
# -----------------------------------------------------------------------


class _QuantumIterativeScheme(qc._IterativeScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(file_name: str) -> _QuantumIterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the initial guess.
        """
        hdf_data = qu.HDF().read(
            file_name,
            [
                qu.HDF.EntryKeys.WVG.value,
                qu.HDF.EntryKeys.SSF.value,
                qu.HDF.EntryKeys.ADR.value,
                qu.HDF.EntryKeys.MATSUBARA.value,
            ],
        )
        return _QuantumIterativeScheme.Guess(
            hdf_data[qu.HDF.EntryKeys.WVG.value],
            hdf_data[qu.HDF.EntryKeys.SSF.value],
            np.ascontiguousarray(hdf_data[qu.HDF.EntryKeys.ADR.value]),
            hdf_data[qu.HDF.EntryKeys.MATSUBARA.value],
        )

    # Save results to disk
    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        if scheme.inputs.degeneracy > 0:
            pd.DataFrame(scheme.adr).to_hdf(
                self.hdf_file_name, key=qu.HDF.EntryKeys.ADR.value
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

        def to_native(self) -> qn.QStlsGuess:
            native_guess = qn.QstlsGuess()
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
        scheme = qn.Qstls(inputs.to_native())
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

        def to_native(self) -> qn.QstlsInput:
            native_input = qn.QstlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.to_native())
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
        self._unpack_fixed_adr_files(inputs)
        scheme = qn.Qstls(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)
        self._zip_fixed_adr_files(inputs)
        self._clean_fixed_adr_files(scheme.inputs)

    # Unpack zip folder with fixed component of the auxiliary density response
    @qu.MPI.run_only_on_root
    def _unpack_fixed_adr_files(self, inputs) -> None:
        fixed_iet_source_file = inputs.fixed_iet
        if inputs.fixed_iet != "":
            inputs.fixed_iet = "qupled_tmp_run_directory"
        if fixed_iet_source_file != "":
            with zipfile.ZipFile(fixed_iet_source_file, "r") as zip_file:
                zip_file.extractall(inputs.fixed_iet)

    # Save results to disk
    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        super()._save(scheme)
        pd.DataFrame(scheme.bf).to_hdf(self.hdf_file_name, key="bf")

    # Zip all files for the fixed component of the auxiliary density response
    @qu.MPI.run_only_on_root
    def _zip_fixed_adr_files(self, inputs) -> None:
        if inputs.fixed_iet == "":
            degeneracy = inputs.degeneracy
            matsubara = inputs.matsubara
            theory = inputs.theory
            adr_file = f"adr_fixed_theta{degeneracy:5.3f}_matsubara{matsubara}_{theory}"
            adr_file_zip = f"{adr_file}.zip"
            adr_file_bin = f"{adr_file}_wv*.bin"
            with zipfile.ZipFile(adr_file_zip, "w") as zip_file:
                for bin_file in glob.glob(adr_file_bin):
                    zip_file.write(bin_file)
                    os.remove(bin_file)

    # Remove temporary run directory
    @qu.MPI.run_only_on_root
    def _clean_fixed_adr_files(self, inputs) -> None:
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

        def to_native(self) -> qn.QstlsInput:
            native_input = qn.QstlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.to_native())
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
        self._unpack_fixed_adr_files(inputs)
        scheme = qn.QVSStls(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)
        self._zip_fixed_adr_files(inputs)
        self._clean_fixed_adr_files(scheme.inputs)

    # Unpack zip folder with fixed component of the auxiliary density response
    @qu.MPI.run_only_on_root
    def _unpack_fixed_adr_files(self, inputs) -> None:
        fixed_source_file = inputs.fixed
        if inputs.fixed != "":
            inputs.fixed = "qupled_tmp_run_directory"
        if fixed_source_file != "":
            with zipfile.ZipFile(fixed_source_file, "r") as zip_file:
                zip_file.extractall(inputs.fixed)

    # Save results to disk
    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        super()._save(scheme)
        pd.DataFrame(scheme.free_energy_grid).to_hdf(
            self.hdf_file_name, key=qu.HDF.EntryKeys.FXC_GRID.value
        )
        pd.DataFrame(scheme.free_energy_integrand).to_hdf(
            self.hdf_file_name, key=qu.HDF.EntryKeys.FXCI.value
        )
        pd.DataFrame(scheme.alpha).to_hdf(
            self.hdf_file_name, key=qu.HDF.EntryKeys.ALPHA.value
        )

    # Zip all files for the fixed component of the auxiliary density response
    @qu.MPI.run_only_on_root
    def _zip_fixed_adr_files(self, inputs) -> None:
        if inputs.fixed == "":
            degeneracy = inputs.degeneracy
            matsubara = inputs.matsubara
            theory = inputs.theory
            adr_file_zip = (
                f"adr_fixed_theta{degeneracy:5.3f}_matsubara{matsubara}_{theory}.zip"
            )
            adr_file_bin = "THETA*.bin"
            with zipfile.ZipFile(adr_file_zip, "w") as zip_file:
                for bin_file in glob.glob(adr_file_bin):
                    zip_file.write(bin_file)
                    os.remove(bin_file)

    # Remove the temporary run directory
    @qu.MPI.run_only_on_root
    def _clean_fixed_adr_files(self, inputs) -> None:
        if os.path.isdir(inputs.fixed):
            shutil.rmtree(inputs.fixed)

    # Set the free energy integrand from a dataframe produced in output
    @staticmethod
    def get_free_energy_integrand(file_name: str) -> qn.FreeEnergyIntegrand:
        return qc.VSStls.get_free_energy_integrand(file_name)

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

        def to_native(self) -> qn.QVSStlsInput:
            native_input = qn.QVSStlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.to_native())
                else:
                    setattr(native_input, attr, value)
            return native_input
