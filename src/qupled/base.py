from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

from . import native
from . import util


class ClassicScheme:

    def __init__(self):
        # File to store output on disk
        self.hdf_file_name: str = None  #: Name of the output file.

    # Compute the scheme
    def _compute(self, scheme) -> None:
        self.hdf_file_name = self._get_hdf_file(scheme.inputs)
        status = scheme.compute()
        self._check_status_and_clean(status, scheme.recovery)

    # Check that the dielectric scheme was solved without errors
    @util.MPI.run_only_on_root
    def _check_status_and_clean(self, status: bool, recovery: str) -> None:
        """Checks that the scheme was solved correctly and removes temporary files generated at run-time

        Args:
            status: status obtained from the native code. If status == 0 the scheme was solved correctly.
            recovery: name of the recovery file.
        """
        if status == 0:
            if os.path.isfile(recovery):
                os.remove(recovery)
            print("Dielectric theory solved successfully!")
        else:
            sys.exit("Error while solving the dielectric theory")

    # Save results to disk
    def _get_hdf_file(self, inputs) -> str:
        """Sets the name of the hdf file used to store the output

        Args:
            inputs: input parameters
        """
        coupling = inputs.coupling
        degeneracy = inputs.degeneracy
        theory = inputs.theory
        return f"rs{coupling:5.3f}_theta{degeneracy:5.3f}_{theory}.h5"

    @util.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        inputs = scheme.inputs
        """Stores the results obtained by solving the scheme."""
        pd.DataFrame(
            {
                util.HDF.EntryKeys.COUPLING.value: inputs.coupling,
                util.HDF.EntryKeys.DEGENERACY.value: inputs.degeneracy,
                util.HDF.EntryKeys.THEORY.value: inputs.theory,
                util.HDF.EntryKeys.RESOLUTION.value: inputs.resolution,
                util.HDF.EntryKeys.CUTOFF.value: inputs.cutoff,
                util.HDF.EntryKeys.FREQUENCY_CUTOFF.value: inputs.frequency_cutoff,
                util.HDF.EntryKeys.MATSUBARA.value: inputs.matsubara,
            },
            index=[util.HDF.EntryKeys.INFO.value],
        ).to_hdf(self.hdf_file_name, key=util.HDF.EntryKeys.INFO.value, mode="w")
        if inputs.degeneracy > 0:
            pd.DataFrame(scheme.idr).to_hdf(
                self.hdf_file_name, key=util.HDF.EntryKeys.IDR.value
            )
            pd.DataFrame(scheme.sdr).to_hdf(
                self.hdf_file_name, key=util.HDF.EntryKeys.SDR.value
            )
            pd.DataFrame(scheme.slfc).to_hdf(
                self.hdf_file_name, key=util.HDF.EntryKeys.SLFC.value
            )
        pd.DataFrame(scheme.ssf).to_hdf(
            self.hdf_file_name, key=util.HDF.EntryKeys.SSF.value
        )
        pd.DataFrame(scheme.ssf_HF).to_hdf(
            self.hdf_file_name, key=util.HDF.EntryKeys.SSF_HF.value
        )
        pd.DataFrame(scheme.wvg).to_hdf(
            self.hdf_file_name, key=util.HDF.EntryKeys.WVG.value
        )

    # Compute radial distribution function
    def compute_rdf(
        self, rdf_grid: np.ndarray = None, write_to_hdf: bool = True
    ) -> np.array:
        """Computes the radial distribution function from the data stored in the output file.

        Args:
            rdf_grid: The grid used to compute the radial distribution function.
                Default = ``None`` (see :func:`qupled.util.Hdf.computeRdf`)
            write_to_hdf: Flag marking whether the rdf data should be added to the output hdf file, default to True

        Returns:
            The radial distribution function

        """
        if util.MPI().rank() > 0:
            write_to_hdf = False
        return util.HDF().compute_rdf(self.hdf_file_name, rdf_grid, write_to_hdf)

    # Compute the internal energy
    def compute_internal_energy(self) -> float:
        """Computes the internal energy from the data stored in the output file.

        Returns:
            The internal energy

        """
        return util.HDF().compute_internal_energy(self.hdf_file_name)

    # Plot results
    @util.MPI.run_only_on_root
    def plot(
        self,
        to_plot: list[str],
        matsubara: np.ndarray = None,
        rdf_grid: np.ndarray = None,
    ) -> None:
        """Plots the results stored in the output file.

        Args:
            to_plot: A list of quantities to plot. This list can include all the values written to the
                 output hdf file. The radial distribution function is computed and added to the output
                 file if necessary
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Default = None, all matsubara frequencies are plotted)
            rdf_grid: The grid used to compute the radial distribution function. Applies only when the radial
                distribution function is plotted. Default = ``None`` (see :func:`qupled.util.Hdf.computeRdf`).

        """
        if util.HDF.EntryKeys.RDF.value in to_plot:
            self.compute_rdf(rdf_grid)
        util.HDF().plot(self.hdf_file_name, to_plot, matsubara)


# -----------------------------------------------------------------------
# IterativeScheme class
# -----------------------------------------------------------------------


class IterativeScheme(ClassicScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(file_name: str) -> IterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the initial guess.
        """
        hdf_data = util.HDF().read(
            file_name, [util.HDF.EntryKeys.WVG.value, util.HDF.EntryKeys.SLFC.value]
        )
        return IterativeScheme.Guess(
            hdf_data[util.HDF.EntryKeys.WVG.value],
            hdf_data[util.HDF.EntryKeys.SLFC.value],
        )

    # Save results to disk
    @util.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        inputs = scheme.inputs
        pd.DataFrame(
            {
                util.HDF.EntryKeys.COUPLING.value: inputs.coupling,
                util.HDF.EntryKeys.DEGENERACY.value: inputs.degeneracy,
                util.HDF.EntryKeys.ERROR.value: scheme.error,
                util.HDF.EntryKeys.THEORY.value: inputs.theory,
                util.HDF.EntryKeys.RESOLUTION.value: inputs.resolution,
                util.HDF.EntryKeys.CUTOFF.value: inputs.cutoff,
                util.HDF.EntryKeys.FREQUENCY_CUTOFF.value: inputs.frequency_cutoff,
                util.HDF.EntryKeys.MATSUBARA.value: inputs.matsubara,
            },
            index=[util.HDF.EntryKeys.INFO.value],
        ).to_hdf(self.hdf_file_name, key=util.HDF.EntryKeys.INFO.value)

    # Initial guess
    class Guess:

        def __init__(self, wvg: np.ndarray = None, slfc: np.ndarray = None):
            self.wvg = wvg
            """ Wave-vector grid. Default = ``None``"""
            self.slfc = slfc
            """ Static local field correction. Default = ``None``"""

        def to_native(self) -> native.StlsGuess:
            native_guess = native.StlsGuess()
            for attr, value in self.__dict__.items():
                native_value = value if value is not None else np.empty(0)
                setattr(native_guess, attr, native_value)
            return native_guess


# -----------------------------------------------------------------------
# QuantumIterativeScheme class
# -----------------------------------------------------------------------


class QuantumIterativeScheme(IterativeScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(file_name: str) -> QuantumIterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the initial guess.
        """
        hdf_data = util.HDF().read(
            file_name,
            [
                util.HDF.EntryKeys.WVG.value,
                util.HDF.EntryKeys.SSF.value,
                util.HDF.EntryKeys.ADR.value,
                util.HDF.EntryKeys.MATSUBARA.value,
            ],
        )
        return QuantumIterativeScheme.Guess(
            hdf_data[util.HDF.EntryKeys.WVG.value],
            hdf_data[util.HDF.EntryKeys.SSF.value],
            np.ascontiguousarray(hdf_data[util.HDF.EntryKeys.ADR.value]),
            hdf_data[util.HDF.EntryKeys.MATSUBARA.value],
        )

    # Save results to disk
    @util.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        if scheme.inputs.degeneracy > 0:
            pd.DataFrame(scheme.adr).to_hdf(
                self.hdf_file_name, key=util.HDF.EntryKeys.ADR.value
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

        def to_native(self) -> native.QStlsGuess:
            native_guess = native.QstlsGuess()
            for attr, value in self.__dict__.items():
                native_value = value if value is not None else np.empty(0)
                setattr(native_guess, attr, native_value)
            return native_guess
