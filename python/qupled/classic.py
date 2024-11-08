from __future__ import annotations
import sys
import os
from abc import ABC, abstractmethod
from typing import Callable
import numpy as np
import pandas as pd
import qupled.util as qu
import qupled.qupled as qp

# -----------------------------------------------------------------------
# ClassicScheme class
# -----------------------------------------------------------------------


class ClassicScheme(ABC):

    # Compute the scheme
    def computeScheme(self, compute: Callable[None, int], save: Callable) -> None:
        status = compute()
        self._checkStatusAndClean(status)
        save()

    # Check that the dielectric scheme was solved without errors
    @qu.MPI.runOnlyOnRoot
    def _checkStatusAndClean(self, status: bool) -> None:
        """Checks that the scheme was solved correctly and removes temporarary files generated at run-time

        Args:
            status: status obtained from the native code. If status == 0 the scheme was solved correctly.
        """
        if status == 0:
            if os.path.isfile(self.recovery):
                os.remove(self.recovery)
            print("Dielectric theory solved successfully!")
        else:
            sys.exit("Error while solving the dielectric theory")

    # Save results to disk
    def _getHdfFile(self) -> str:
        """Sets the name of the hdf file used to store the output"""
        return "rs%5.3f_theta%5.3f_%s.h5" % (
            self.inputs.coupling,
            self.inputs.degeneracy,
            self.inputs.theory,
        )

    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """Stores the results obtained by solving the scheme."""
        pd.DataFrame(
            {
                "coupling": self.inputs.coupling,
                "degeneracy": self.inputs.degeneracy,
                "theory": self.inputs.theory,
                "resolution": self.inputs.resolution,
                "cutoff": self.inputs.cutoff,
                "matsubara": self.inputs.matsubara,
            },
            index=["info"],
        ).to_hdf(self.hdfFileName, key="info", mode="w")
        pd.DataFrame(self.idr).to_hdf(self.hdfFileName, key="idr")
        pd.DataFrame(self.sdr).to_hdf(self.hdfFileName, key="sdr")
        pd.DataFrame(self.slfc).to_hdf(self.hdfFileName, key="slfc")
        pd.DataFrame(self.ssf).to_hdf(self.hdfFileName, key="ssf")
        pd.DataFrame(self.ssfHF).to_hdf(self.hdfFileName, key="ssfHF")
        pd.DataFrame(self.wvg).to_hdf(self.hdfFileName, key="wvg")

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
        return qp.computeInternalEnergy(self.wvg, self.ssf, self.inputs.coupling)

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


class RpaMetaclass(type(ClassicScheme), type(qp.Rpa)):
    pass


class Rpa(qp.Rpa, ClassicScheme, metaclass=RpaMetaclass):
    """
    Args:
        inputs: Input parameters.
    """

    # Constructor
    def __init__(self, inputs: Rpa.Input):
        # Construct the base classes
        super().__init__(inputs)
        # File to store output on disk
        self.hdfFileName: str = self._getHdfFile()  #: Name of the output file.

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """
        Solves the scheme and saves the results.
        """
        super().computeScheme(super().compute, self._save)

    # Input class
    class Input(qp.RpaInput):
        """
        Class used to manage the input for the :obj:`qupled.classic.Rpa` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__()
            self.coupling: float = coupling
            self.degeneracy: float = degeneracy
            self.chemicalPotential: list[float] = [-10.0, 10.0]
            self.matsubara: int = 128
            self.resolution: float = 0.1
            self.cutoff: float = 10.0
            self.intError: float = 1.0e-5
            self.int2DScheme: str = "full"
            self.threads: int = 1
            self.theory: str = "RPA"

        @property
        def coupling(self) -> float:
            """Coupling parameter."""
            return super().coupling

        @property
        def degeneracy(self) -> float:
            """Degeneracy parameter."""
            return super().degeneracy

        @property
        def chemicalPotential(self) -> list[float]:
            """
            Initial guess for the chemical potential. Default = ``[-10, 10]``
            """
            return super().chemicalPotential

        @property
        def matsubara(self) -> int:
            """Number of Matsubara frequencies. Default = ``128``"""
            return super().matsubara

        @property
        def resolution(self) -> float:
            """Resolution of the wave-vector grid. Default = ``0.1``"""
            return super().resolution

        @property
        def cutoff(self) -> float:
            """Cutoff for the wave-vector grid. Default = ``10.0``"""
            return super().cutoff

        @property
        def intError(self) -> float:
            """Accuracy (relative error) in the computation of integrals. Default = ``1.0e-5``"""
            return super().intError

        @property
        def int2DScheme(self) -> str:
            """
            Scheme used to solve two-dimensional integrals
            allowed options include:

            - full: the inner integral is evaluated at arbitrary points
              selected automatically by the quadrature rule

            - segregated: the inner integral is evaluated on a fixed
              grid that depends on the integrand that is being processed

            Segregated is usually faster than full but it could become
            less accurate if the fixed points are not chosen correctly. Default = ``'full'``
            """
            return super().int2DScheme

        @property
        def threads(self) -> int:
            """Number of OMP threads for parallel calculations. Default = ``1``"""
            return super().threads

        # Setters

        @coupling.setter
        def coupling(self, value: float):
            super(Rpa.Input, self.__class__).coupling.fset(self, value)

        @degeneracy.setter
        def degeneracy(self, value: float):
            super(Rpa.Input, self.__class__).degeneracy.fset(self, value)

        @chemicalPotential.setter
        def chemicalPotential(self, value: list[float]):
            super(Rpa.Input, self.__class__).chemicalPotential.fset(self, value)

        @matsubara.setter
        def matsubara(self, value: int):
            super(Rpa.Input, self.__class__).matsubara.fset(self, value)

        @resolution.setter
        def resolution(self, value: float):
            super(Rpa.Input, self.__class__).resolution.fset(self, value)

        @cutoff.setter
        def cutoff(self, value: float):
            super(Rpa.Input, self.__class__).cutoff.fset(self, value)

        @intError.setter
        def intError(self, value: float):
            super(Rpa.Input, self.__class__).intError.fset(self, value)

        @int2DScheme.setter
        def int2DScheme(self, value: str):
            super(Rpa.Input, self.__class__).int2DScheme.fset(self, value)

        @threads.setter
        def threads(self, value: int):
            super(Rpa.Input, self.__class__).threads.fset(self, value)


# -----------------------------------------------------------------------
# ESA class
# -----------------------------------------------------------------------


class ESAMetaclass(type(ClassicScheme), type(qp.ESA)):
    pass


class ESA(ClassicScheme, qp.ESA, metaclass=ESAMetaclass):
    """
    Args:
        inputs: Input parameters.
    """

    # Constructor
    def __init__(self, inputs: ESA.Input):
        # Construct the base classes
        super().__init__(inputs)
        # File to store output on disk
        self.hdfFileName: str = self._getHdfFile()  #: Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """
        Solves the scheme and saves the results.
        """
        super().computeScheme(super().compute, self._save)

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
# IterativeScheme class
# -----------------------------------------------------------------------


class IterativeScheme(ClassicScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def getInitialGuess(fileName: str) -> qp.StlsGuess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            fileName : name of the file used to extract the information for the initial guess.
        """
        guess = qp.StlsGuess()
        hdfData = qu.Hdf().read(fileName, ["wvg", "slfc"])
        guess.wvg = hdfData["wvg"]
        guess.slfc = hdfData["slfc"]
        return guess

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save()
        pd.DataFrame(
            {
                "coupling": self.inputs.coupling,
                "degeneracy": self.inputs.degeneracy,
                "error": self.error,
                "theory": self.inputs.theory,
                "resolution": self.inputs.resolution,
                "cutoff": self.inputs.cutoff,
                "matsubara": self.inputs.matsubara,
            },
            index=["info"],
        ).to_hdf(self.hdfFileName, key="info")


# -----------------------------------------------------------------------
# Stls class
# -----------------------------------------------------------------------


class StlsMetaclass(type(IterativeScheme), type(qp.Stls)):
    pass


class Stls(IterativeScheme, qp.Stls, metaclass=StlsMetaclass):
    """
    Args:
        inputs: Input parameters.
    """

    # Constructor
    def __init__(self, inputs: Stls.Input):
        # Construct the base classes
        super().__init__(inputs)
        # File to store output on disk
        self.hdfFileName: str = self._getHdfFile()  #: Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """
        Solves the scheme and saves the results.
        """
        super().computeScheme(super().compute, self._save)

    # Input class
    class Input(Rpa.Input, qp.StlsInput):
        """
        Class used to manage the input for the :obj:`qupled.classic.Stls` class.
        """

        def __init__(self, coupling: float, degeneracy: float, initGuess: bool = True):
            super().__init__(coupling, degeneracy)
            self.error: float = 1.0e-5
            self.mixing: float = 1.0
            self.iterations: int = 1000
            self.outputFrequency: int = 10
            self.recoveryFile: str = ""
            if initGuess:
                self.guess: qp.StlsGuess = qp.StlsGuess()
            self.theory: str = "STLS"

        @property
        def error(self) -> float:
            """Minimum error for convergence. Default = ``1.0e-5``"""
            return super().error

        @property
        def mixing(self) -> float:
            """Mixing parameter. Default = ``1.0``"""
            return super().mixing

        @property
        def iterations(self) -> int:
            """Maximum number of iterations. Default = ``1000``"""
            return super().iterations

        @property
        def outputFrequency(self) -> int:
            """Output frequency to write the recovery file. Default = ``10``"""
            return super().outputFrequency

        @property
        def recoveryFile(self) -> str:
            """Name of the recovery file. Default = ``""``"""
            return super().recoveryFile

        @property
        def guess(self) -> qupled.qupled.StlsGuess:
            """Initial guess."""
            return super().guess

        # Setters

        @error.setter
        def error(self, value: float):
            super(Stls.Input, self.__class__).error.fset(self, value)

        @mixing.setter
        def mixing(self, value: float):
            super(Stls.Input, self.__class__).mixing.fset(self, value)

        @iterations.setter
        def iterations(self, value: int):
            super(Stls.Input, self.__class__).iterations.fset(self, value)

        @outputFrequency.setter
        def outputFrequency(self, value: int):
            super(Stls.Input, self.__class__).outputFrequency.fset(self, value)

        @recoveryFile.setter
        def recoveryFile(self, value: str):
            super(Stls.Input, self.__class__).recoveryFile.fset(self, value)

        @guess.setter
        def guess(self, value: qupled.qupled.StlsGuess):
            super(Stls.Input, self.__class__).guess.fset(self, value)


# -----------------------------------------------------------------------
# StlsIet class
# -----------------------------------------------------------------------


class StlsIet(IterativeScheme, qp.Stls, metaclass=StlsMetaclass):
    """
    Args:
        inputs: Input parameters.
    """

    # Constructor
    def __init__(self, inputs: StlsIet.Input):
        # Construct the base classes
        super().__init__(inputs)
        # File to store output on disk
        self.hdfFileName: str = self._getHdfFile()  #: Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results."""
        super().computeScheme(super().compute, self._save)

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save()
        pd.DataFrame(self.bf).to_hdf(self.hdfFileName, key="bf")

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

        @property
        def mapping(self) -> str:
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
            return super().mapping

        @mapping.setter
        def mapping(self, value: str):
            super(StlsIet.Input, self.__class__).mapping.fset(self, value)


# -----------------------------------------------------------------------
# VSStls class
# -----------------------------------------------------------------------


class VSStlsMetaclass(type(IterativeScheme), type(qp.VSStls)):
    pass


class VSStls(IterativeScheme, qp.VSStls, metaclass=VSStlsMetaclass):
    """
    Args:
        inputs: Input parameters.
    """

    # Constructor
    def __init__(self, inputs: VSStls.Input):
        # Construct the base classes
        super().__init__(inputs)
        # File to store output on disk
        self.hdfFileName: str = self._getHdfFile()  #: Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results."""
        super().computeScheme(super().compute, self._save)

    # Save results
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save()
        pd.DataFrame(self.freeEnergyGrid).to_hdf(self.hdfFileName, key="fxcGrid")
        pd.DataFrame(self.freeEnergyIntegrand).to_hdf(self.hdfFileName, key="fxci")
        pd.DataFrame(self.alpha).to_hdf(self.hdfFileName, key="alpha")

    # Set the free energy integrand from a dataframe produced in output
    @staticmethod
    def getFreeEnergyIntegrand(fileName: str) -> qp.FreeEnergyIntegrand():
        """Constructs the free energy integrand by extracting the information from an output file.

        Args:
            fileName : name of the file used to extract the information for the free energy integrand.
        """
        fxci = qp.FreeEnergyIntegrand()
        hdfData = qu.Hdf().read(fileName, ["fxcGrid", "fxci", "alpha"])
        fxci.grid = hdfData["fxcGrid"]
        fxci.integrand = np.ascontiguousarray(hdfData["fxci"])
        fxci.alpha = hdfData["alpha"]
        return fxci

    # Input class
    class Input(Stls.Input, qp.VSStlsInput):
        """
        Class used to manage the input for the :obj:`qupled.classic.VSStls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.alpha: list[float] = [0.5, 1.0]
            self.couplingResolution: float = 0.1
            self.degeneracyResolution: float = 0.1
            self.errorAlpha: float = 1.0e-3
            self.iterationsAlpha: int = 50
            self.freeEnergyIntegrand: qupled.FreeEnergyIntegrand = (
                qp.FreeEnergyIntegrand()
            )
            self.threads: int = 9
            self.theory: str = "VSSTLS"

        @property
        def alpha(self) -> list[float]:
            """Initial guess for the free parameter. Default = ``[0.5, 1.0]``"""
            return super().alpha

        @property
        def couplingResolution(self) -> float:
            """Resolution of the coupling parameter grid. Default = ``0.1``"""
            return super().couplingResolution

        @property
        def degeneracyResolution(self) -> float:
            """Resolution of the degeneracy parameter grid. Default = ``0.1``"""
            return super().degeneracyResolution

        @property
        def errorAlpha(self) -> float:
            """Minimum error for convergence in the free parameter. Default = ``1.0e-3``"""
            return super().errorAlpha

        @property
        def iterationsAlpha(self) -> int:
            """Maximum number of iterations to determine the free parameter. Default = ``50``"""
            return super().iterationsAlpha

        @property
        def freeEnergyIntegrand(self) -> qupled.FreeEnergyIntegrand:
            """Pre-computed free energy integrand."""
            return super().freeEnergyIntegrand

        @property
        def threads(self) -> int:
            """Number of threads. Default = ``9``"""
            return super().threads

        # Setters

        @alpha.setter
        def alpha(self, value: list[float]):
            super(VSStls.Input, self.__class__).alpha.fset(self, value)

        @couplingResolution.setter
        def couplingResolution(self, value: float):
            super(VSStls.Input, self.__class__).couplingResolution.fset(self, value)

        @degeneracyResolution.setter
        def degeneracyResolution(self, value: float):
            super(VSStls.Input, self.__class__).degeneracyResolution.fset(self, value)

        @errorAlpha.setter
        def errorAlpha(self, value: float):
            super(VSStls.Input, self.__class__).errorAlpha.fset(self, value)

        @iterationsAlpha.setter
        def iterationsAlpha(self, value: int):
            super(VSStls.Input, self.__class__).iterationsAlpha.fset(self, value)

        @freeEnergyIntegrand.setter
        def freeEnergyIntegrand(self, value: qupled.FreeEnergyIntegrand):
            super(VSStls.Input, self.__class__).freeEnergyIntegrand.fset(self, value)

        @threads.setter
        def threads(self, value: int):
            super(VSStls.Input, self.__class__).threads.fset(self, value)
