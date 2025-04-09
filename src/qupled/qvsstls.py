from __future__ import annotations

import glob
import os
import shutil
import zipfile

from . import mpi
from . import native
from . import qstls
from . import vsstls


class QVSStls(vsstls.VSStls):
    """
    Class used to solve the QVStls scheme.
    """

    def __init__(self):
        super().__init__()
        self.results: Result = Result()
        # Undocumented properties
        self.native_scheme_cls = native.QVSStls
        self.native_inputs = native.QVSStlsInput()

    # Compute
    def compute(self, inputs: Input):
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self._unpack_fixed_adr_files(inputs)
        super().compute(inputs)
        self._zip_fixed_adr_files(inputs)
        self._clean_fixed_adr_files(inputs)

    @staticmethod
    def get_initial_guess(run_id: int, database_name: str | None = None) -> qstls.Guess:
        return qstls.Qstls.get_initial_guess(run_id, database_name)

    # Unpack zip folder with fixed component of the auxiliary density response
    @mpi.MPI.run_only_on_root
    def _unpack_fixed_adr_files(self, inputs):
        """
        Unpacks a fixed adr file into a temporary directory.

        This method extracts the contents of a fixed source file (if provided)
        into a temporary directory named "qupled_tmp_run_directory". If the
        `inputs.fixed` attribute is not an empty string, it is updated to point
        to the temporary directory.

        Args:
            inputs: An object containing the `fixed` attribute, which represents
                    the path to the fixed source file. If `inputs.fixed` is an
                    empty string, no action is taken.

        Raises:
            zipfile.BadZipFile: If the provided fixed source file is not a valid
                                ZIP file.
            FileNotFoundError: If the fixed source file does not exist.
        """
        fixed_source_file = inputs.fixed
        if inputs.fixed != "":
            inputs.fixed = "qupled_tmp_run_directory"
        if fixed_source_file != "":
            with zipfile.ZipFile(fixed_source_file, "r") as zip_file:
                zip_file.extractall(inputs.fixed)

    # Zip all files for the fixed component of the auxiliary density response
    @mpi.MPI.run_only_on_root
    def _zip_fixed_adr_files(self, inputs):
        """
        Compresses and removes binary files matching a specific pattern into a ZIP archive.

        This method creates a ZIP file containing binary files with names matching the
        pattern "THETA*.bin". The name of the ZIP file is generated based on the
        `degeneracy`, `matsubara`, and `theory` attributes of the `inputs` object.
        After adding the binary files to the ZIP archive, the original binary files
        are deleted from the filesystem.

        Args:
            inputs: An object containing the following attributes:
                - fixed (str): A string that determines if the operation should proceed.
                  If empty, the method executes; otherwise, it does nothing.
                - degeneracy (float): A value used to format the ZIP file name.
                - matsubara (int): A value used to format the ZIP file name.
                - theory (str): A string used to format the ZIP file name.
        """
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
    @mpi.MPI.run_only_on_root
    def _clean_fixed_adr_files(self, inputs):
        """
        Removes the directory specified by the `fixed` attribute of the `inputs` object.

        This method checks if the path provided in `inputs.fixed` is a directory.
        If it is, the directory and all its contents are deleted.

        Args:
            inputs: An object that contains a `fixed` attribute, which is the path
                    to the directory to be removed.

        Raises:
            OSError: If the directory cannot be removed due to permission issues
                     or if the path is invalid.
        """
        if os.path.isdir(inputs.fixed):
            shutil.rmtree(inputs.fixed)


# Input class
class Input(vsstls.Input, qstls.Input):
    """
    Class used to manage the input for the :obj:`qupled.qvsstls.QVSStls` class.
    """

    def __init__(self, coupling: float, degeneracy: float):
        vsstls.Input.__init__(self, coupling, degeneracy)
        qstls.Input.__init__(self, coupling, degeneracy)
        # Undocumented default values
        self.theory: str = "QVSSTLS"


# Result
class Result(vsstls.Result, qstls.Result):
    """
    Class used to manage the results for the :obj:`qupled.qvsstls.QVSStls` class.
    """

    def __init__(self):
        super().__init__()
