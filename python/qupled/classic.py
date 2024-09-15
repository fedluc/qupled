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

    # Setup inputs object
    def _setInputs(
        self,
        coupling: float,
        degeneracy: float,
        theory: str,
        chemicalPotential: list[float],
        cutoff: float,
        matsubara: int,
        resolution: float,
    ) -> None:
        # """ Sets up the content of :obj:`inputs` """
        self.inputs.coupling = coupling
        self.inputs.degeneracy = degeneracy
        self.inputs.theory = theory
        self.inputs.chemicalPotential = chemicalPotential
        self.inputs.cutoff = cutoff
        self.inputs.matsubara = matsubara
        self.inputs.resolution = resolution
        self.inputs.intError = 1.0e-5
        self.inputs.threads = 1

    # Compute the scheme
    @abstractmethod
    def compute(self):
        pass

    # Check input before computing
    def _checkInputs(self) -> None:
        """Checks that the content of :obj:`inputs` is correct"""
        if self.inputs.theory not in self.allowedTheories:
            sys.exit("Invalid dielectric theory")

    # Check that the dielectric scheme was solved without errors
    @qu.MPI.runOnlyOnRoot
    def _checkStatusAndClean(self, status: bool) -> None:
        """Checks that the scheme was solved correctly and removes temporarary files generated at run-time

        Args:
            status: status obtained from the native code. If status == 0 the scheme was solved correctly.
        """
        if status == 0:
            if os.path.isfile(self.scheme.recovery):
                os.remove(self.scheme.recovery)
            print("Dielectric theory solved successfully!")
        else:
            sys.exit("Error while solving the dielectric theory")

    # Save results to disk
    def _setHdfFile(self) -> None:
        """Sets the name of the hdf file used to store the output"""
        self.hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (
            self.inputs.coupling,
            self.inputs.degeneracy,
            self.inputs.theory,
        )

    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """Stores the results obtained by solving the scheme."""
        assert self.scheme is not None
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
        pd.DataFrame(self.scheme.idr).to_hdf(self.hdfFileName, key="idr")
        pd.DataFrame(self.scheme.sdr).to_hdf(self.hdfFileName, key="sdr")
        pd.DataFrame(self.scheme.slfc).to_hdf(self.hdfFileName, key="slfc")
        pd.DataFrame(self.scheme.ssf).to_hdf(self.hdfFileName, key="ssf")
        pd.DataFrame(self.scheme.ssfHF).to_hdf(self.hdfFileName, key="ssfHF")
        pd.DataFrame(self.scheme.wvg).to_hdf(self.hdfFileName, key="wvg")

    # Compute radial distribution function
    def computeRdf(
        self, rdfGrid: np.ndarray = None, writeToHdf: bool = True
    ) -> np.array:
        """Computes the radial distribution function from the data stored in the output file.

        Args:
            rdfGrid: The grid used to compute the radial distribution function.
                (Defaults to None, see :func:`qupled.util.Hdf.computeRdf`)
            writeToHdf: Flag marking whether the rdf data should be added to the output hdf file, default to True

        Returns:
            The radial distribution function

        """
        self._checkSolution("compute the radial distribution function")
        if qu.MPI().getRank() > 0:
            writeToHdf = False
        return qu.Hdf().computeRdf(self.hdfFileName, rdfGrid, writeToHdf)

    # Compute the internal energy
    def computeInternalEnergy(self) -> float:
        """Computes the internal energy from the data stored in the output file.

        Returns:
            The internal energy

        """
        self._checkSolution("compute the internal energy")
        return qp.computeInternalEnergy(
            self.scheme.wvg, self.scheme.ssf, self.inputs.coupling
        )

    # Plot results
    @qu.MPI.runOnlyOnRoot
    def plot(
        self,
        toPlot: list[str],
        matsubara: np.ndarray = None,
        rdfGrid: np.ndarray = None,
    ) -> None:
        """Plots the results stored in the output file`.

        Args:
            toPlot: A list of quantities to plot. This list can include all the values written to the
                 output hdf file. The radial distribution funciton is computed and added to the output
                 file if necessary
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Default = None, all matsubara frequencies are plotted)
            rdfGrid: The grid used to compute the radial distribution function. Applies only when the radial
                distribution function is plotted (Default = None, see :func:`qupled.classic.Stls.computeRdf`)

        """
        self._checkSolution("plot results")
        if "rdf" in toPlot:
            self.computeRdf(rdfGrid)
        qu.Hdf().plot(self.hdfFileName, toPlot, matsubara)

    # Check if a solution is available to perform a given action
    def _checkSolution(self, action: str) -> None:
        """Check if a solution is available to be used

        Args:
            action: Name of the action to be performed. Only used to print an error message if no solution is found

        """
        if self.scheme is None:
            sys.exit("No solution to " + action)


# -----------------------------------------------------------------------
# ClassicSchemeNew class
# -----------------------------------------------------------------------


class ClassicSchemeNew(ABC):

    # Compute the scheme
    def compute(self, compute: Callable[None, int], save: Callable) -> None:
        self._checkInputs()
        status = compute()
        self._checkStatusAndClean(status)
        self._setHdfFile()
        save()

    # Check input before computing
    def _checkInputs(self) -> None:
        """Checks that the content of :obj:`inputs` is correct"""
        if self.inputs.theory not in self.allowedTheories:
            sys.exit("Invalid dielectric theory")

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
    def _setHdfFile(self) -> None:
        """Sets the name of the hdf file used to store the output"""
        self.hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (
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
                (Defaults to None, see :func:`qupled.util.Hdf.computeRdf`)
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
        """Plots the results stored in the output file`.

        Args:
            toPlot: A list of quantities to plot. This list can include all the values written to the
                 output hdf file. The radial distribution funciton is computed and added to the output
                 file if necessary
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Default = None, all matsubara frequencies are plotted)
            rdfGrid: The grid used to compute the radial distribution function. Applies only when the radial
                distribution function is plotted (Default = None, see :func:`qupled.classic.Stls.computeRdf`)

        """
        if "rdf" in toPlot:
            self.computeRdf(rdfGrid)
        qu.Hdf().plot(self.hdfFileName, toPlot, matsubara)


# -----------------------------------------------------------------------
# RPA class
# -----------------------------------------------------------------------


class RpaMetaclass(type(ClassicSchemeNew), type(qp.Rpa)):
    pass


class Rpa(qp.Rpa, ClassicSchemeNew, metaclass=RpaMetaclass):
    """
    Class used to setup and solve the classical Randon-Phase approximaton scheme as described by
    `Bohm and Pines <https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.609>`_.
    The inputs used to solve the scheme are defined when creating the class, but can be
    later modified by changing the attribute :obj:`inputs`. After the solution is completed
    the results are saved to an hdf file and can be plotted via the method :obj:`plot`.

    Args:
        inputs: Input parameters.
    """

    # Constructor
    def __init__(self, inputs: qp.RpaInput):
        # Construct the base classes
        super().__init__(inputs)
        # Allowed theories
        self.allowedTheories = ["RPA"]
        # Input object
        self.inputs: qp.RpaInput = inputs
        # File to store output on disk
        self.hdfFileName: str = None  #: Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results.

        The results are stored as pandas dataframes in an hdf file with the following keywords:

        - info: A dataframe containing information on the input parameters, it includes:

          - coupling: the coupling parameter,
          - degeneracy: the degeneracy parameter,
          - theory: the theory that is being solved,
          - resolution: the resolution in the wave-vector grid,
          - cutoff: the cutoff in the wave-vector grid,
          - matsubara: the number of matsubara frequencies

        - idr (*ndarray*, 2D): the ideal density response
        - sdr (*ndarray*):  the static density response
        - slfc (*ndarray*):  the static local field correction
        - ssf (*ndarray*):  the static structure factor
        - ssfHF (*ndarray*):  the Hartree-Fock static structure factor
        - wvg (*ndarray*):  the wave-vector grid

        If the radial distribution function was computed (see computeRdf), then the hdf file contains
        two additional keywords:

        - rdf (*ndarray*):  the radial distribution function
        - rdfGrid (*ndarray*):  the grid used to compute the radial distribution function

        The name of the hdf file is stored in :obj:`hdfFileName`.
        """
        ClassicSchemeNew.compute(self, qp.Rpa.compute, self._save)


# -----------------------------------------------------------------------
# ESA class
# -----------------------------------------------------------------------


class ESAMetaclass(type(ClassicSchemeNew), type(qp.Rpa)):
    pass


class ESA(ClassicSchemeNew, qp.ESA, metaclass=ESAMetaclass):
    """
    Class used to setup and solve the Effective Static Approximation scheme as described by
    `Dornheim and collaborators <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.165102>`_.
    The inputs used to solve the scheme are defined when creating the class, but can be
    later modified by changing the attribute :obj:`inputs`. After the solution is completed
    the results are saved to an hdf file and can be plotted via the method :obj:`plot`.

    Args:
        inputs: Input parameters used to solve the scheme.
    """

    # Constructor
    def __init__(self, inputs: qp.RpaInput):
        # Construct the base classes
        super().__init__(inputs)
        # Allowed theories
        self.allowedTheories = ["ESA"]
        # Input object
        self.inputs: qp.RpaInput = inputs
        # File to store output on disk
        self.hdfFileName: str = None  #: Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results.

        The results are stored as pandas dataframes in an hdf file with the following keywords:

        - info: A dataframe containing information on the input parameters, it includes:

          - coupling: the coupling parameter,
          - degeneracy: the degeneracy parameter,
          - theory: the theory that is being solved,
          - resolution: the resolution in the wave-vector grid,
          - cutoff: the cutoff in the wave-vector grid,
          - matsubara: the number of matsubara frequencies

        - idr (*ndarray*, 2D): the ideal density response
        - sdr (*ndarray*):  the static density response
        - slfc (*ndarray*):  the static local field correction
        - ssf (*ndarray*):  the static structure factor
        - ssfHF (*ndarray*):  the Hartree-Fock static structure factor
        - wvg (*ndarray*):  the wave-vector grid

        If the radial distribution function was computed (see computeRdf), then the hdf file contains
        two additional keywords:

        - rdf (*ndarray*):  the radial distribution function
        - rdfGrid (*ndarray*):  the grid used to compute the radial distribution function

        The name of the hdf file is stored in :obj:`hdfFileName`.
        """
        ClassicSchemeNew.compute(self, qp.ESA.compute, self._save)


# -----------------------------------------------------------------------
# Stls class
# -----------------------------------------------------------------------


class StlsMetaclass(type(ClassicSchemeNew), type(qp.Stls)):
    pass


class Stls(ClassicSchemeNew, qp.Stls, metaclass=StlsMetaclass):
    """
    Class used to setup and solve the classical STLS scheme as described by
    `Tanaka and Ichimaru <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_.
    The inputs used to solve the scheme are defined when creating the class, but can be
    later modified by changing the attribute :obj:`inputs`. After the solution is completed
    the results are saved to an hdf file and can be plotted via the method :obj:`plot`.

    Args:
        inputs: Input parameters used to solve the scheme.
    """

    # Constructor
    def __init__(self, inputs: qp.StlsInput):
        # Construct the base classes
        super().__init__(inputs)
        # Allowed theories
        self.allowedTheories = ["STLS"]
        # Input object
        self.inputs: qp.StlsInput = inputs
        # File to store output on disk
        self.hdfFileName: str = None  #: Name of the output file

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results.

        The results are stored as pandas dataframes in an hdf file with the following keywords:

        - info: A dataframe containing information on the input parameters, it includes:

          - coupling: the coupling parameter,
          - degeneracy: the degeneracy parameter,
          - error: the residual error at the end of the solution
          - theory: the theory that is being solved,
          - resolution: the resolution in the wave-vector grid,
          - cutoff: the cutoff in the wave-vector grid,
          - matsubara: the number of matsubara frequencies

        - idr (*ndarray*, 2D): the ideal density response
        - sdr (*ndarray*):  the static density response
        - slfc (*ndarray*):  the static local field correction
        - ssf (*ndarray*):  the static structure factor
        - ssfHF (*ndarray*):  the Hartree-Fock static structure factor
        - wvg (*ndarray*):  the wave-vector grid

        If the radial distribution function was computed (see computeRdf), then the hdf file contains
        two additional keywords:

        - rdf (*ndarray*):  the radial distribution function
        - rdfGrid (*ndarray*):  the grid used to compute the radial distribution function

        The name of the hdf file is stored in :obj:`hdfFileName`.
        """
        ClassicSchemeNew.compute(self, qp.Stls.compute, self._save)

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

    # Set the initial guess from a dataframe produced in output
    def getInitialGuess(self, fileName: str) -> qp.StlsGuess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            fileName : name of the file used to extract the information for the initial guess.
        """
        guess = qp.StlsGuess()
        hdfData = qu.Hdf().read(fileName, ["wvg", "slfc"])
        guess.wvg = hdfData["wvg"]
        guess.slfc = hdfData["slfc"]
        return guess


# -----------------------------------------------------------------------
# StlsIet class
# -----------------------------------------------------------------------


class StlsIet(Stls):
    """

    Class used to setup and solve the classical STLS-IET scheme as described by
    `Tanaka <https://aip.scitation.org/doi/full/10.1063/1.4969071>`_ and by
    `Tolias and collaborators <https://aip.scitation.org/doi/full/10.1063/1.4969071>`_.
    This class inherits most of its methods and attributes from :obj:`qupled.classic.Stls`.

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.
        chemicalPotential: Initial guess for the chemical potential, defaults to [-10.0, 10.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for convergence, defaults to 1.0e-5.
        mapping: Classical to quantum mapping. See :func:`qupled.qupled.StlsInput.iet`
        mixing: Mixing parameter for iterative solution, defaults to 1.0.
        guess:  Initial guess for the iterative solution, defaults to None, i.e. slfc = 0.
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
        scheme2DIntegrals: numerical scheme used to solve two-dimensional integrals. See :func:`qupled.qupled.Input.int2DScheme`
    """

    # Constructor
    def __init__(
        self,
        coupling: float,
        degeneracy: float,
        theory: str,
        chemicalPotential: list[float] = [-10.0, 10.0],
        cutoff: float = 10.0,
        error: float = 1.0e-5,
        mapping: str = "standard",
        mixing: float = 1.0,
        guess: qp.StlsGuess = None,
        iterations: int = 1000,
        matsubara: int = 128,
        outputFrequency: int = 10,
        recoveryFile: str = None,
        resolution: float = 0.1,
        scheme2DIntegrals: str = "full",
    ):
        # Allowed theories
        self.allowedTheories = ["STLS-HNC", "STLS-IOI", "STLS-LCT"]
        # Input object
        self.inputs: qupled.qupled.StlsInput = (
            qp.StlsInput()
        )  #: Inputs to solve the scheme.
        super()._setInputs(
            coupling,
            degeneracy,
            theory,
            chemicalPotential,
            cutoff,
            error,
            mixing,
            guess,
            iterations,
            matsubara,
            outputFrequency,
            recoveryFile,
            resolution,
        )
        self.inputs.iet = mapping
        self.inputs.int2DScheme = scheme2DIntegrals
        self._checkInputs()
        # Scheme to solve
        self.scheme: qp.Stls = None
        # File to store output on disk
        self.hdfFileName = None

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results to and hdf file. Extends the output produced by
        :func:`qupled.classic.Stls.compute` by adding the option to save the bridge function adder
        as a new dataframe in the hdf file. The bridge function adder dataframe can be accessed as `bf`.
        """
        self._checkInputs()
        self.scheme = qp.Stls(self.inputs)
        status = self.scheme.compute()
        self._checkStatusAndClean(status)
        self._setHdfFile()
        self._save()

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save()
        pd.DataFrame(self.scheme.bf).to_hdf(self.hdfFileName, key="bf")


# -----------------------------------------------------------------------
# VSStls class
# -----------------------------------------------------------------------


class VSStls(Stls):
    """
    Class used to setup and solve the classical VS-STLS scheme as described by
    `Vashishta and Singwi <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.6.875>`_ and by
    `Sjostrom and Dufty <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.115123>`_.
    This class inherits most of its methods and attributes from :obj:`qupled.classic.Stls`

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.
        chemicalPotential: Initial guess for the chemical potential, defaults to [-100.0, 100.0].
        cutoff:  Cutoff for the wave-vector grid, defaults to 10.0.
        error: Minimum error for convergence, defaults to 1.0e-5.
        mixing: Mixing parameter for iterative solution, defaults to 1.0.
        iterations: Maximum number of iterations, defaults to 1000.
        matsubara: Number of matsubara frequencies, defaults to 128.
        outputFrequency: Frequency used to print the recovery files, defaults to 10.
        recoveryFile: Name of the recovery file used to restart the simulation, defualts to None.
        resolution: Resolution of the wave-vector grid, defaults to 0.1.
        alpha: Initial guess for the free parameter, defaults to [0.5, 1.0]
        couplingResolution: Resolution of the coupling parameter grid, defaults to 0.01
        degeneracyResolution: Resolution of the degeneracy parameter grid, defaults to 0.01
        errorAlpha: Minimum error for convergence in the free parameter iterations, defaults to 1.0e-3
        iterationsAlpha: Maximum number of iterations for the free parameter, defaults to 50
        errorIntegrals: Accuracy (as a relative error) for the integral computations, defaults to 1.0-5
        threads: number of OMP threads for parallel calculations, defualts to 1
    """

    # Constructor
    def __init__(
        self,
        coupling: float,
        degeneracy: float,
        chemicalPotential: list[float] = [-100.0, 100.0],
        cutoff: float = 10.0,
        error: float = 1.0e-5,
        mixing: float = 1.0,
        guess: qp.StlsGuess = None,
        iterations: int = 1000,
        matsubara: int = 128,
        outputFrequency: int = 10,
        recoveryFile: str = None,
        resolution: float = 0.1,
        alpha: list[float] = [0.5, 1.0],
        couplingResolution: float = 0.01,
        degeneracyResolution: float = 0.01,
        errorAlpha: float = 1.0e-3,
        iterationsAlpha: int = 50,
        errorIntegrals: float = 1.0e-5,
        threads: int = 1,
    ):
        # Allowed theories
        self.allowedTheories: list[str] = ["VSSTLS"]
        # Input object
        self.inputs: qupled.qupled.VSStlsInput = (
            qp.VSStlsInput()
        )  #: Inputs to solve the scheme.
        super()._setInputs(
            coupling,
            degeneracy,
            "VSSTLS",
            chemicalPotential,
            cutoff,
            error,
            mixing,
            guess,
            iterations,
            matsubara,
            outputFrequency,
            recoveryFile,
            resolution,
        )
        self.inputs.alpha = alpha
        self.inputs.couplingResolution = couplingResolution
        self.inputs.degeneracyResolution = degeneracyResolution
        self.inputs.errorAlpha = errorAlpha
        self.inputs.iterationsAlpha = iterationsAlpha
        self.inputs.threads = threads
        self.inputs.intError = errorIntegrals
        # Scheme to solve
        self.scheme: qp.VSStls = None
        # File to store output on disk
        self.hdfFileName = None

    # Compute
    @qu.MPI.recordTime
    @qu.MPI.synchronizeRanks
    def compute(self) -> None:
        """Solves the scheme and saves the results to and hdf file. Extends the output produced by
        :func:`qupled.classic.Stls.compute` by adding the option to save the free energy integrand
        and the corresponding coupling parameter grid as a new dataframe in the hdf file. The free
        energy integrand dataframe can be accessed as `fxci` and the corresponding coupling parameter
        grid data frame as `fxcGrid`.
        """
        self._checkInputs()
        self.scheme = qp.VSStls(self.inputs)
        status = self.scheme.compute()
        self._checkStatusAndClean(status)
        self._setHdfFile()
        self._save()

    # Save results
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save()
        pd.DataFrame(self.scheme.freeEnergyGrid).to_hdf(self.hdfFileName, key="fxcGrid")
        pd.DataFrame(self.scheme.freeEnergyIntegrand).to_hdf(
            self.hdfFileName, key="fxci"
        )
        pd.DataFrame(self.scheme.alpha).to_hdf(self.hdfFileName, key="alpha")

    # Set the free energy integrand from a dataframe produced in output
    def setFreeEnergyIntegrand(self, fileName: str) -> None:
        """Constructs the free energy integrand by extracting the information from an output file.

        Args:
            fileName : name of the file used to extract the information for the free energy integrand.
        """
        fxci = qp.FreeEnergyIntegrand()
        hdfData = qu.Hdf().read(fileName, ["fxcGrid", "fxci", "alpha"])
        fxci.grid = hdfData["fxcGrid"]
        fxci.integrand = np.ascontiguousarray(hdfData["fxci"])
        fxci.alpha = hdfData["alpha"]
        self.inputs.freeEnergyIntegrand = fxci
