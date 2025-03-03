from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

from qupled import native
import qupled.util as qu


# -----------------------------------------------------------------------
# _ClassicScheme class
# -----------------------------------------------------------------------


class _ClassicScheme:

    def __init__(self):
        # File to store output on disk
        self.hdf_file_name: str = None  #: Name of the output file.

    # Compute the scheme
    def _compute(self, scheme) -> None:
        self.hdf_file_name = self._get_hdf_file(scheme.inputs)
        status = scheme.compute()
        self._check_status_and_clean(status, scheme.recovery)

    # Check that the dielectric scheme was solved without errors
    @qu.MPI.run_only_on_root
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

    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        inputs = scheme.inputs
        """Stores the results obtained by solving the scheme."""
        pd.DataFrame(
            {
                qu.HDF.EntryKeys.COUPLING.value: inputs.coupling,
                qu.HDF.EntryKeys.DEGENERACY.value: inputs.degeneracy,
                qu.HDF.EntryKeys.THEORY.value: inputs.theory,
                qu.HDF.EntryKeys.RESOLUTION.value: inputs.resolution,
                qu.HDF.EntryKeys.CUTOFF.value: inputs.cutoff,
                qu.HDF.EntryKeys.FREQUENCY_CUTOFF.value: inputs.frequency_cutoff,
                qu.HDF.EntryKeys.MATSUBARA.value: inputs.matsubara,
            },
            index=[qu.HDF.EntryKeys.INFO.value],
        ).to_hdf(self.hdf_file_name, key=qu.HDF.EntryKeys.INFO.value, mode="w")
        if inputs.degeneracy > 0:
            pd.DataFrame(scheme.idr).to_hdf(
                self.hdf_file_name, key=qu.HDF.EntryKeys.IDR.value
            )
            pd.DataFrame(scheme.sdr).to_hdf(
                self.hdf_file_name, key=qu.HDF.EntryKeys.SDR.value
            )
            pd.DataFrame(scheme.slfc).to_hdf(
                self.hdf_file_name, key=qu.HDF.EntryKeys.SLFC.value
            )
        pd.DataFrame(scheme.ssf).to_hdf(
            self.hdf_file_name, key=qu.HDF.EntryKeys.SSF.value
        )
        pd.DataFrame(scheme.ssf_HF).to_hdf(
            self.hdf_file_name, key=qu.HDF.EntryKeys.SSF_HF.value
        )
        pd.DataFrame(scheme.wvg).to_hdf(
            self.hdf_file_name, key=qu.HDF.EntryKeys.WVG.value
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
        if qu.MPI().rank() > 0:
            write_to_hdf = False
        return qu.HDF().compute_rdf(self.hdf_file_name, rdf_grid, write_to_hdf)

    # Compute the internal energy
    def compute_internal_energy(self) -> float:
        """Computes the internal energy from the data stored in the output file.

        Returns:
            The internal energy

        """
        return qu.HDF().compute_internal_energy(self.hdf_file_name)

    # Plot results
    @qu.MPI.run_only_on_root
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
        if qu.HDF.EntryKeys.RDF.value in to_plot:
            self.compute_rdf(rdf_grid)
        qu.HDF().plot(self.hdf_file_name, to_plot, matsubara)


# -----------------------------------------------------------------------
# RPA class
# -----------------------------------------------------------------------


class Rpa(_ClassicScheme):

    # Compute
    @qu.MPI.record_time
    @qu.MPI.synchronize_ranks
    def compute(self, inputs: Rpa.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.Rpa(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)

    # Input class
    class Input:
        """
        Class used to manage the input for the :obj:`qupled.classic.Rpa` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            self.coupling: float = coupling
            """Coupling parameter."""
            self.degeneracy: float = degeneracy
            """Degeneracy parameter."""
            self.chemical_potential: list[float] = [-10.0, 10.0]
            """Initial guess for the chemical potential. Default = ``[-10, 10]``"""
            self.matsubara: int = 128
            """Number of Matsubara frequencies. Default = ``128``"""
            self.resolution: float = 0.1
            """Resolution of the wave-vector grid. Default =  ``0.1``"""
            self.cutoff: float = 10.0
            """Cutoff for the wave-vector grid. Default =  ``10.0``"""
            self.frequency_cutoff: float = 10.0
            """Cutoff for the frequency (applies only in the ground state). Default =  ``10.0``"""
            self.integral_error: float = 1.0e-5
            """Accuracy (relative error) in the computation of integrals. Default = ``1.0e-5``"""
            self.integral_strategy: str = "full"
            """
            Scheme used to solve two-dimensional integrals
            allowed options include:

            - full: the inner integral is evaluated at arbitrary points
              selected automatically by the quadrature rule

            - segregated: the inner integral is evaluated on a fixed
              grid that depends on the integrand that is being processed

            Segregated is usually faster than full but it could become
            less accurate if the fixed points are not chosen correctly. Default =  ``'full'``
            """
            self.threads: int = 1
            """Number of OMP threads for parallel calculations. Default =  ``1``"""
            self.theory: str = "RPA"

        def to_native(self) -> native.RpaInput:
            native_input = native.RpaInput()
            for attr, value in self.__dict__.items():
                setattr(native_input, attr, value)
            return native_input


# -----------------------------------------------------------------------
# ESA class
# -----------------------------------------------------------------------


class ESA(_ClassicScheme):
    """
    Args:
        inputs: Input parameters.
    """

    # Compute
    @qu.MPI.record_time
    @qu.MPI.synchronize_ranks
    def compute(self, inputs: ESA.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.ESA(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)

    # Input class
    class Input(Rpa.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.ESA` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            # Undocumented default values
            self.theory = "ESA"


# -----------------------------------------------------------------------
# _IterativeScheme class
# -----------------------------------------------------------------------


class _IterativeScheme(_ClassicScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(file_name: str) -> _IterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the initial guess.
        """
        hdf_data = qu.HDF().read(
            file_name, [qu.HDF.EntryKeys.WVG.value, qu.HDF.EntryKeys.SLFC.value]
        )
        return _IterativeScheme.Guess(
            hdf_data[qu.HDF.EntryKeys.WVG.value],
            hdf_data[qu.HDF.EntryKeys.SLFC.value],
        )

    # Save results to disk
    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        inputs = scheme.inputs
        pd.DataFrame(
            {
                qu.HDF.EntryKeys.COUPLING.value: inputs.coupling,
                qu.HDF.EntryKeys.DEGENERACY.value: inputs.degeneracy,
                qu.HDF.EntryKeys.ERROR.value: scheme.error,
                qu.HDF.EntryKeys.THEORY.value: inputs.theory,
                qu.HDF.EntryKeys.RESOLUTION.value: inputs.resolution,
                qu.HDF.EntryKeys.CUTOFF.value: inputs.cutoff,
                qu.HDF.EntryKeys.FREQUENCY_CUTOFF.value: inputs.frequency_cutoff,
                qu.HDF.EntryKeys.MATSUBARA.value: inputs.matsubara,
            },
            index=[qu.HDF.EntryKeys.INFO.value],
        ).to_hdf(self.hdf_file_name, key=qu.HDF.EntryKeys.INFO.value)

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
# Stls class
# -----------------------------------------------------------------------


class Stls(_IterativeScheme):

    # Compute
    @qu.MPI.record_time
    @qu.MPI.synchronize_ranks
    def compute(self, inputs: Stls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.Stls(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)

    # Input class
    class Input(Rpa.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.Stls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.error: float = 1.0e-5
            """Minimum error for convergence. Default = ``1.0e-5``"""
            self.mixing: float = 1.0
            """Mixing parameter. Default = ``1.0``"""
            self.iterations: int = 1000
            """Maximum number of iterations. Default = ``1000``"""
            self.output_frequency: int = 10
            """Output frequency to write the recovery file. Default = ``10``"""
            self.recovery_file: str = ""
            """Name of the recovery file. Default = ``""``"""
            self.guess: Stls.Guess = Stls.Guess()
            """Initial guess. Default = ``Stls.Guess()``"""
            # Undocumented default values
            self.theory: str = "STLS"

        def to_native(self) -> native.StlsInput:
            native_input = native.StlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.to_native())
                else:
                    setattr(native_input, attr, value)
            return native_input


# -----------------------------------------------------------------------
# StlsIet class
# -----------------------------------------------------------------------


class StlsIet(_IterativeScheme):

    # Compute
    @qu.MPI.record_time
    @qu.MPI.synchronize_ranks
    def compute(self, inputs: StlsIet.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.Stls(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)

    # Save results to disk
    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        pd.DataFrame(scheme.bf).to_hdf(self.hdf_file_name, key="bf")

    # Input class
    class Input(Stls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.StlsIet` class.
        Accepted theories: ``STLS-HNC``, ``STLS-IOI`` and ``STLS-LCT``.
        """

        def __init__(self, coupling: float, degeneracy: float, theory: str):
            super().__init__(coupling, degeneracy)
            if theory not in {"STLS-HNC", "STLS-IOI", "STLS-LCT"}:
                sys.exit("Invalid dielectric theory")
            self.theory = theory
            self.mapping = "standard"
            r"""
            Mapping for the classical-to-quantum coupling parameter
            :math:`\Gamma` used in the iet schemes. Allowed options include:

            - standard: :math:`\Gamma \propto \Theta^{-1}`

            - sqrt: :math:`\Gamma \propto (1 + \Theta)^{-1/2}`

            - linear: :math:`\Gamma \propto (1 + \Theta)^{-1}`

            where :math:`\Theta` is the degeneracy parameter. Far from the ground state
            (i.e. :math:`\Theta\gg1`) all mappings lead identical results, but at
            the ground state they can differ significantly (the standard
            mapping diverges). Default = ``standard``.
            """

        def to_native(self) -> native.StlsInput:
            native_input = native.StlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.to_native())
                else:
                    setattr(native_input, attr, value)
            return native_input


# -----------------------------------------------------------------------
# VSStls class
# -----------------------------------------------------------------------


class VSStls(_IterativeScheme):

    # Compute
    @qu.MPI.record_time
    @qu.MPI.synchronize_ranks
    def compute(self, inputs: VSStls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.VSStls(inputs.to_native())
        self._compute(scheme)
        self._save(scheme)

    # Save results
    @qu.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
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

    # Set the free energy integrand from a dataframe produced in output
    @staticmethod
    def get_free_energy_integrand(file_name: str) -> native.FreeEnergyIntegrand:
        """Constructs the free energy integrand by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the free energy integrand.
        """
        fxci = native.FreeEnergyIntegrand()
        hdf_data = qu.HDF().read(
            file_name,
            [
                qu.HDF.EntryKeys.FXC_GRID.value,
                qu.HDF.EntryKeys.FXCI.value,
                qu.HDF.EntryKeys.ALPHA.value,
            ],
        )
        fxci.grid = hdf_data[qu.HDF.EntryKeys.FXC_GRID.value]
        fxci.integrand = np.ascontiguousarray(hdf_data[qu.HDF.EntryKeys.FXCI.value])
        fxci.alpha = hdf_data[qu.HDF.EntryKeys.ALPHA.value]
        return fxci

    # Input class
    class Input(Stls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.VSStls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.alpha: list[float] = [0.5, 1.0]
            """Initial guess for the free parameter. Default = ``[0.5, 1.0]``"""
            self.coupling_resolution: float = 0.1
            """Resolution of the coupling parameter grid. Default = ``0.1``"""
            self.degeneracy_resolution: float = 0.1
            """Resolution of the degeneracy parameter grid. Default = ``0.1``"""
            self.error_alpha: float = 1.0e-3
            """Minimum error for convergence in the free parameter. Default = ``1.0e-3``"""
            self.iterations_alpha: int = 50
            """Maximum number of iterations to determine the free parameter. Default = ``50``"""
            self.free_energy_integrand: native.FreeEnergyIntegrand = (
                native.FreeEnergyIntegrand()
            )
            """Pre-computed free energy integrand."""
            self.threads: int = 9
            """Number of threads. Default = ``9``"""
            # Undocumented default values
            self.theory: str = "VSSTLS"

        def to_native(self) -> native.VSStlsInput:
            native_input = native.VSStlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.to_native())
                else:
                    setattr(native_input, attr, value)
            return native_input
