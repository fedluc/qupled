# -----------------------------------------------------------------------
# QstlsIet class
# -----------------------------------------------------------------------

from __future__ import annotations
import glob
import os
import shutil
import zipfile

from . import native
from . import util
from . import base
from . import qstls
from . import stlsiet


class QstlsIet(base.QuantumIterativeScheme):
    """
    Args:
        inputs: Input parameters.
    """

    # Compute
    def compute(self, inputs: QstlsIet.Input):
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self._unpack_fixed_adr_files(inputs)
        super().compute(inputs, native.Qstls, native.QstlsInput(), self.Result())
        self._zip_fixed_adr_files(inputs)
        self._clean_fixed_adr_files(inputs)

    # Unpack zip folder with fixed component of the auxiliary density response
    @util.MPI.run_only_on_root
    def _unpack_fixed_adr_files(self, inputs):
        """
        Unpacks a fixed IET source file into a temporary directory.

        This method extracts the contents of a fixed IET source file, if provided, into a temporary
        directory named "qupled_tmp_run_directory". The method updates the `inputs.fixed_iet` attribute
        to point to this temporary directory.

        Args:
            inputs: An object containing the `fixed_iet` attribute, which is the path to the fixed IET
                    source file. If the attribute is an empty string, no action is taken.
        """
        fixed_iet_source_file = inputs.fixed_iet
        if inputs.fixed_iet != "":
            inputs.fixed_iet = "qupled_tmp_run_directory"
        if fixed_iet_source_file != "":
            with zipfile.ZipFile(fixed_iet_source_file, "r") as zip_file:
                zip_file.extractall(inputs.fixed_iet)

    # Zip all files for the fixed component of the auxiliary density response
    @util.MPI.run_only_on_root
    def _zip_fixed_adr_files(self, inputs):
        """
        Compresses and removes binary files matching a specific naming pattern into a ZIP archive.

        This method generates a filename based on the provided `inputs` object, which contains
        attributes such as `degeneracy`, `matsubara`, and `theory`. If the `fixed_iet` attribute
        of `inputs` is an empty string, it constructs a filename pattern for binary files and
        compresses all matching files into a ZIP archive. After adding the files to the archive,
        the original binary files are deleted.

        Args:
            inputs: An object containing the following attributes:
                - fixed_iet (str): A string that determines whether the operation proceeds.
                - degeneracy (float): A numerical value used in the filename.
                - matsubara (int): An integer value used in the filename.
                - theory (str): A string value used in the filename.

        Raises:
            FileNotFoundError: If no files matching the binary file pattern are found.
            OSError: If there is an issue writing to the ZIP file or removing the binary files.
        """
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
    @util.MPI.run_only_on_root
    def _clean_fixed_adr_files(self, inputs):
        """
        Removes the directory specified by `inputs.fixed_iet` if it exists.

        This method checks if the path provided in `inputs.fixed_iet` is a directory.
        If it is, the directory and all its contents are deleted.

        Args:
            inputs: An object containing the attribute `fixed_iet`, which is the
                    path to the directory to be removed.
        """
        if os.path.isdir(inputs.fixed_iet):
            shutil.rmtree(inputs.fixed_iet)

    # Input class
    class Input(stlsiet.StlsIet.Input, qstls.Qstls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.QStlsIet` class.
        Accepted theories: ``QSTLS-HNC``, ``QSTLS-IOI`` and ``QSTLS-LCT``.
        """

        def __init__(self, coupling: float, degeneracy: float, theory: str):
            stlsiet.StlsIet.Input.__init__(self, coupling, degeneracy, "STLS-HNC")
            qstls.Qstls.Input.__init__(self, coupling, degeneracy)
            if theory not in {"QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"}:
                raise ValueError("Invalid dielectric theory")
            self.theory = theory
            self.fixed_iet = ""
            """
            Name of the zip file storing the iet part of the fixed components
            of the auxiliary density response. Default = ``""``
            """

    # Result class
    class Result(stlsiet.StlsIet.Result, qstls.Qstls.Result):
        """
        Class used to store the results for the :obj:`qupled.quantum.QstlsIet` class.
        """

        def __init__(self):
            super().__init__()
