# -----------------------------------------------------------------------
# QVSStls class
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
from . import vsstls


class QVSStls(base.QuantumIterativeScheme):

    # Compute
    def compute(self, inputs: QVSStls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        self._unpack_fixed_adr_files(inputs)
        super().compute(inputs, native.QVSStls, native.QVSStlsInput(), self.Result())
        self._zip_fixed_adr_files(inputs)
        self._clean_fixed_adr_files(inputs)

    # Unpack zip folder with fixed component of the auxiliary density response
    @util.MPI.run_only_on_root
    def _unpack_fixed_adr_files(self, inputs) -> None:
        fixed_source_file = inputs.fixed
        if inputs.fixed != "":
            inputs.fixed = "qupled_tmp_run_directory"
        if fixed_source_file != "":
            with zipfile.ZipFile(fixed_source_file, "r") as zip_file:
                zip_file.extractall(inputs.fixed)

    # Zip all files for the fixed component of the auxiliary density response
    @util.MPI.run_only_on_root
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
    @util.MPI.run_only_on_root
    def _clean_fixed_adr_files(self, inputs) -> None:
        if os.path.isdir(inputs.fixed):
            shutil.rmtree(inputs.fixed)

    # Set the free energy integrand from a dataframe produced in output
    @staticmethod
    def get_free_energy_integrand(file_name: str) -> native.FreeEnergyIntegrand:
        return vsstls.VSStls.get_free_energy_integrand(file_name)

    # Input class
    class Input(vsstls.VSStls.Input, qstls.Qstls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.QVSStls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            vsstls.VSStls.Input.__init__(self, coupling, degeneracy)
            qstls.Qstls.Input.__init__(self, coupling, degeneracy)
            # Undocumented default values
            self.theory: str = "QVSSTLS"

    # Result
    class Result(vsstls.VSStls.Result, qstls.Qstls.Result):

        def __init__(self):
            super().__init__()
