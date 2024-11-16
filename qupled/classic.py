from __future__ import annotations
import sys
import os
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
        self.hdfFileName: str = None  #: Name of the output file.

    # Compute the scheme
    def _compute(self, scheme) -> None:
        self.hdfFileName = self._getHdfFile(scheme.inputs)
        status = scheme.compute()
        self._checkStatusAndClean(status, scheme.recovery)

    # Check that the dielectric scheme was solved without errors
    @qu.MPI.runOnlyOnRoot
    def _checkStatusAndClean(self, status: bool, recovery: str) -> None:
        """Checks that the scheme was solved correctly and removes temporarary files generated at run-time

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
    def _getHdfFile(self, inputs) -> str:
        """Sets the name of the hdf file used to store the output

        Args:
            inputs: input parameters
        """
        coupling = inputs.coupling
        degeneracy = inputs.degeneracy
        theory = inputs.theory
        return f"rs{coupling:5.3f}_theta{degeneracy:5.3f}_{theory}.h5"

    @qu.MPI.runOnlyOnRoot
    def _save(self, scheme) -> None:
        inputs = scheme.inputs
        """Stores the results obtained by solving the scheme."""
        pd.DataFrame(
            {
                "coupling": inputs.coupling,
                "degeneracy": inputs.degeneracy,
                "theory": inputs.theory,
                "resolution": inputs.resolution,
                "cutoff": inputs.cutoff,
                "matsubara": inputs.matsubara,
            },
            index=["info"],
        ).to_hdf(self.hdfFileName, key="info", mode="w")
        pd.DataFrame(scheme.idr).to_hdf(self.hdfFileName, key="idr")
        pd.DataFrame(scheme.sdr).to_hdf(self.hdfFileName, key="sdr")
        pd.DataFrame(scheme.slfc).to_hdf(self.hdfFileName, key="slfc")
        pd.DataFrame(scheme.ssf).to_hdf(self.hdfFileName, key="ssf")
        pd.DataFrame(scheme.ssfHF).to_hdf(self.hdfFileName, key="ssfHF")
        pd.DataFrame(scheme.wvg).to_hdf(self.hdfFileName, key="wvg")

    # Compute radial distribution function
    def computeRdf(
        self, rdfGrid: np.ndarray = None, writeToHdf: bool = True
    ) -> np.array:
        """Computes the radial distribution function from the data stored in the output file.

        Args:
            rdfGrid: The grid used to compute the radial distribution function.
                Default = ``None`` (see :func:`qupled.util.Hdf.computeRdf`)
            writeToHdf: Flag marking whether the rdf data should be added to the output hdf file, default to True

        Returns:
            The radial distribution function

        """
        if qu.MPI().getRank() > 0:
            writeToHdf = False
        return qu.Hdf().computeRdf(self.hdfFileName, rdfGrid, writeToHdf)

    # Compute the internal energy
    def computeInternalEnergy(self) -> float:
        """Computes the internal energy from the data stored in the output file.

        Returns:
            The internal energy

        """
        return qu.Hdf().computeInternalEnergy(self.hdfFileName)

    # Plot results
    @qu.MPI.runOnlyOnRoot
    def plot(
        self,
        toPlot: list[str],
        matsubara: np.ndarray = None,
        rdfGrid: np.ndarray = None,
    ) -> None:
        """Plots the results stored in the output file.

        Args:
            toPlot: A list of quantities to plot. This list can include all the values written to the
                 output hdf file. The radial distribution funciton is computed and added to the output
                 file if necessary
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Default = None, all matsubara frequencies are plotted)
            rdfGrid: The grid used to compute the radial distribution function. Applies only when the radial
                distribution function is plotted. Default = ``None`` (see :func:`qupled.util.Hdf.computeRdf`).

        """
        if "rdf" in toPlot:
            self.computeRdf(rdfGrid)
        qu.Hdf().plot(self.hdfFileName, toPlot, matsubara)


# -----------------------------------------------------------------------
# RPA class
# -----------------------------------------------------------------------


class Rpa(_ClassicScheme):

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self, inputs: Rpa.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.Rpa(inputs.toNative())
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
            self.chemicalPotential: list[float] = [-10.0, 10.0]
            """Initial guess for the chemical potential. Default = ``[-10, 10]``"""
            self.matsubara: int = 128
            """Number of Matsubara frequencies. Default = ``128``"""
            self.resolution: float = 0.1
            """Resolution of the wave-vector grid. Default =  ``0.1``"""
            self.cutoff: float = 10.0
            """Cutoff for the wave-vector grid. Default =  ``10.0``"""
            self.intError: float = 1.0e-5
            """Accuracy (relative error) in the computation of integrals. Default = ``1.0e-5``"""
            self.int2DScheme: str = "full"
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

        def toNative(self) -> native.RpaInput:
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
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self, inputs: ESA.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.ESA(inputs.toNative())
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
    def getInitialGuess(fileName: str) -> _IterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            fileName : name of the file used to extract the information for the initial guess.
        """
        hdfData = qu.Hdf().read(fileName, ["wvg", "slfc"])
        return _IterativeScheme.Guess(hdfData["wvg"], hdfData["slfc"])

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        inputs = scheme.inputs
        pd.DataFrame(
            {
                "coupling": inputs.coupling,
                "degeneracy": inputs.degeneracy,
                "error": scheme.error,
                "theory": inputs.theory,
                "resolution": inputs.resolution,
                "cutoff": inputs.cutoff,
                "matsubara": inputs.matsubara,
            },
            index=["info"],
        ).to_hdf(self.hdfFileName, key="info")

    # Initial guess
    class Guess:

        def __init__(self, wvg: np.ndarray = None, slfc: np.ndarray = None):
            self.wvg = wvg
            """ Wave-vector grid. Default = ``None``"""
            self.slfc = slfc
            """ Static local field correction. Default = ``None``"""

        def toNative(self) -> native.StlsGuess:
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
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self, inputs: Stls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.Stls(inputs.toNative())
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
            self.outputFrequency: int = 10
            """Output frequency to write the recovery file. Default = ``10``"""
            self.recoveryFile: str = ""
            """Name of the recovery file. Default = ``""``"""
            self.guess: Stls.Guess = Stls.Guess()
            """Initial guess. Default = ``Stls.Guess()``"""
            # Undocumented default values
            self.theory: str = "STLS"

        def toNative(self) -> native.StlsInput:
            native_input = native.StlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.toNative())
                else:
                    setattr(native_input, attr, value)
            return native_input


# -----------------------------------------------------------------------
# StlsIet class
# -----------------------------------------------------------------------


class StlsIet(_IterativeScheme):

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self, inputs: StlsIet.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.Stls(inputs.toNative())
        self._compute(scheme)
        self._save(scheme)

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        pd.DataFrame(scheme.bf).to_hdf(self.hdfFileName, key="bf")

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

        def toNative(self) -> native.StlsInput:
            native_input = native.StlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.toNative())
                else:
                    setattr(native_input, attr, value)
            return native_input


# -----------------------------------------------------------------------
# VSStls class
# -----------------------------------------------------------------------


class VSStls(_IterativeScheme):

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self, inputs: VSStls.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.VSStls(inputs.toNative())
        self._compute(scheme)
        self._save(scheme)

    # Save results
    @qu.MPI.runOnlyOnRoot
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        pd.DataFrame(scheme.freeEnergyGrid).to_hdf(self.hdfFileName, key="fxcGrid")
        pd.DataFrame(scheme.freeEnergyIntegrand).to_hdf(self.hdfFileName, key="fxci")
        pd.DataFrame(scheme.alpha).to_hdf(self.hdfFileName, key="alpha")

    # Set the free energy integrand from a dataframe produced in output
    @staticmethod
    def getFreeEnergyIntegrand(fileName: str) -> native.FreeEnergyIntegrand():
        """Constructs the free energy integrand by extracting the information from an output file.

        Args:
            fileName : name of the file used to extract the information for the free energy integrand.
        """
        fxci = native.FreeEnergyIntegrand()
        hdfData = qu.Hdf().read(fileName, ["fxcGrid", "fxci", "alpha"])
        fxci.grid = hdfData["fxcGrid"]
        fxci.integrand = np.ascontiguousarray(hdfData["fxci"])
        fxci.alpha = hdfData["alpha"]
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
            self.couplingResolution: float = 0.1
            """Resolution of the coupling parameter grid. Default = ``0.1``"""
            self.degeneracyResolution: float = 0.1
            """Resolution of the degeneracy parameter grid. Default = ``0.1``"""
            self.errorAlpha: float = 1.0e-3
            """Minimum error for convergence in the free parameter. Default = ``1.0e-3``"""
            self.iterationsAlpha: int = 50
            """Maximum number of iterations to determine the free parameter. Default = ``50``"""
            self.freeEnergyIntegrand: qupled.FreeEnergyIntegrand = (
                native.FreeEnergyIntegrand()
            )
            """Pre-computed free energy integrand."""
            self.threads: int = 9
            """Number of threads. Default = ``9``"""
            # Undocumented default values
            self.theory: str = "VSSTLS"

        def toNative(self) -> native.VSStlsInput:
            native_input = native.VSStlsInput()
            for attr, value in self.__dict__.items():
                if attr == "guess":
                    setattr(native_input, attr, value.toNative())
                else:
                    setattr(native_input, attr, value)
            return native_input
