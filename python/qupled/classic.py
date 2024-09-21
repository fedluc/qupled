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


class RpaMetaclass(type(ClassicScheme), type(qp.Rpa)):
    pass


class Rpa(qp.Rpa, ClassicScheme, metaclass=RpaMetaclass):
    """
    Class used to setup and solve the classical Randon-Phase approximaton scheme as described by
    `Bohm and Pines <https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.609>`_.
    The inputs used to solve the scheme are defined when creating the class, but can be
    later modified by changing the attribute :obj:`inputs`. After the solution is completed
    the results are saved to an hdf file and can be plotted via the method :obj:`plot`.

    Args:
        inputs: Input parameters.
    """

    class Input(qp.RpaInput):
        """
        Class used to manage the input for the :obj:`qupled.classic.Rpa` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__()
            self.coupling: float = coupling
            """ Coupling parameter """
            self.degeneracy: float = degeneracy
            """ Degeneracy parameter """
            self.chemicalPotential: list[float] = [-10.0, 10.0]
            """ Initial guess for the chemical potential """
            self.matsubara: int = 128
            """ Number of matsubara frequencies"""
            self.resolution: float = 0.1
            """ Resolution of the wave-vector grid """
            self.cutoff: float = 10.0
            """ cutoff for the wave-vector grid """
            self.intError: float = 1.0e-5
            """ Accuracy (expressed as a relative error) in the computation of the integrals """
            self.int2DScheme: str = "full"
            """ Scheme used to solve two-dimensional integrals
            allowed options include:
        
            - full: the inner integral is evaluated at arbitrary points 
	      selected automatically by the quadrature rule
        
	    - segregated: the inner integral is evaluated on a fixed 
	      grid that depends on the integrand that is being processed
        
            Segregated is usually faster than full but it could become 
	    less accurate if the fixed points are not chosen correctly
            """
            self.threads: int = 1
            """ Number of OMP threads for parallel calculations"""
            # Undocumented default values
            self.theory: list[str] = "RPA"

    # Constructor
    def __init__(self, inputs: Rpa.Input):
        # Construct the base classes
        super().__init__(inputs)
        # File to store output on disk
        self.hdfFileName: str = self._getHdfFile()  #: Name of the output file

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
        super().computeScheme(super().compute, self._save)


# -----------------------------------------------------------------------
# ESA class
# -----------------------------------------------------------------------


class ESAMetaclass(type(ClassicScheme), type(qp.ESA)):
    pass


class ESA(ClassicScheme, qp.ESA, metaclass=ESAMetaclass):
    """
    Class used to setup and solve the Effective Static Approximation scheme as described by
    `Dornheim and collaborators <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.165102>`_.
    The inputs used to solve the scheme are defined when creating the class, but can be
    later modified by changing the attribute :obj:`inputs`. After the solution is completed
    the results are saved to an hdf file and can be plotted via the method :obj:`plot`.

    Args:
        inputs: Input parameters used to solve the scheme.
    """

    class Input(Rpa.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.ESA` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            # Undocumented default values
            self.theory = "ESA"

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
        super().computeScheme(super().compute, self._save)


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
    Class used to setup and solve the classical STLS scheme as described by
    `Tanaka and Ichimaru <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_.
    The inputs used to solve the scheme are defined when creating the class, but can be
    later modified by changing the attribute :obj:`inputs`. After the solution is completed
    the results are saved to an hdf file and can be plotted via the method :obj:`plot`.

    Args:
        inputs: Input parameters used to solve the scheme.
    """

    class Input(Rpa.Input, qp.StlsInput):
        """
        Class used to manage the input for the :obj:`qupled.classic.Stls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            self.error: float = 1.0e-5
            """ minimum error for convergence """
            self.mixing: float = 1.0
            """ mixing paramter """
            self.iterations: int = 1000
            """ Maximum number of iterations """
            self.outputFrequency: int = 10
            """ Output frequency to write the recovery file """
            self.recoveryFile: str = ""
            """ Name of the recovery file """
            self.guess: qp.StlsGuess = qp.StlsGuess()
            """ Initial guess """
            # Undocumented default values
            self.theory = "STLS"

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
        super().computeScheme(super().compute, self._save)


# -----------------------------------------------------------------------
# StlsIet class
# -----------------------------------------------------------------------


class StlsIet(IterativeScheme, qp.Stls, metaclass=StlsMetaclass):
    """

    Class used to setup and solve the classical STLS-IET scheme as described by
    `Tanaka <https://aip.scitation.org/doi/full/10.1063/1.4969071>`_ and by
    `Tolias and collaborators <https://aip.scitation.org/doi/full/10.1063/1.4969071>`_.
    This class inherits most of its methods and attributes from :obj:`qupled.classic.Stls`.

    Args:
        inputs: Input parameters used to solve the scheme.
    """

    class Input(Stls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.StlsIet` class.
        """

        def __init__(self, coupling: float, degeneracy: float, theory: str):
            super().__init__(coupling, degeneracy)
            if theory not in ["STLS-HNC", "STLS-IOI" or "STLS-LCT"]:
                sys.exit("Invalid dielectric theory")
            self.theory = theory
            """ Dielectric theory to solve  """

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
        - bf (*ndarray*): the bridge function adder

        If the radial distribution function was computed (see computeRdf), then the hdf file contains
        two additional keywords:

        - rdf (*ndarray*):  the radial distribution function
        - rdfGrid (*ndarray*):  the grid used to compute the radial distribution function

        The name of the hdf file is stored in :obj:`hdfFileName`.
        """
        super().computeScheme(super().compute, self._save)

    # Save results to disk
    @qu.MPI.runOnlyOnRoot
    def _save(self) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save()
        pd.DataFrame(self.bf).to_hdf(self.hdfFileName, key="bf")


# -----------------------------------------------------------------------
# VSStls class
# -----------------------------------------------------------------------


class VSStlsMetaclass(type(IterativeScheme), type(qp.VSStls)):
    pass


class VSStls(IterativeScheme, qp.VSStls, metaclass=VSStlsMetaclass):
    """
     Class used to setup and solve the classical VS-STLS scheme as described by
    `Vashishta and Singwi <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.6.875>`_ and by
    `Sjostrom and Dufty <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.115123>`_.

    Args:
        inputs: Input parameters used to solve the scheme.
    """

    class Input(Stls.Input, qp.VSStlsInput):
        """
        Class used to manage the input for the :obj:`qupled.classic.Stls` class.
        """

        def __init__(self, coupling: float, degeneracy: float):
            super().__init__(coupling, degeneracy)
            """ Name of the theory that is solved """
            self.alpha: list[float] = [0.5, 1.0]
            """ Initial guess for the free parameter """
            self.couplingResolution: float = 0.1
            """ Resolution of the coupling parameter grid """
            self.degeneracyResolution: float = 0.1
            """ Resolution of the degeneracy parameter grid """
            self.errorAlpha: float = 1.0e-3
            """ Minimum error for convergence in the free parameter """
            self.iterationsAlpha: int = 50
            """ Maximum number of iterations to determine the free parameter """
            self.freeEnergyIntegrand: qupled.FreeEnergyIntegrand = (
                qp.FreeEnergyIntegrand()
            )
            """ Pre-computed free energy integrand """
            # Undocumented default values
            self.threads = 9
            self.theory = "VSSTLS"

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

        - fxcGrid (*ndarray*): coupling parameter grid
        - fxci (*ndarray*): the free energy integrand
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
