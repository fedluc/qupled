# -----------------------------------------------------------------------
# QstlsIet class
# -----------------------------------------------------------------------

from __future__ import annotations
import glob
import os
import shutil
import sys
import zipfile

import pandas as pd

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
    def compute(self, inputs: QstlsIet.Input) -> None:
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
    def _unpack_fixed_adr_files(self, inputs) -> None:
        fixed_iet_source_file = inputs.fixed_iet
        if inputs.fixed_iet != "":
            inputs.fixed_iet = "qupled_tmp_run_directory"
        if fixed_iet_source_file != "":
            with zipfile.ZipFile(fixed_iet_source_file, "r") as zip_file:
                zip_file.extractall(inputs.fixed_iet)

    # Zip all files for the fixed component of the auxiliary density response
    @util.MPI.run_only_on_root
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
    @util.MPI.run_only_on_root
    def _clean_fixed_adr_files(self, inputs) -> None:
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

    class Result(stlsiet.StlsIet.Result, qstls.Qstls.Result):

        def __init__(self):
            super().__init__()
