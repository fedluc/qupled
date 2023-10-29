# Placeholders used to document the c++ classes exposed with python::boost

from __future__ import annotations
import numpy as np

class Input():
    """Base class to handle the inputs.

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.
        theory: theory to solve.
    """
    def __init__(self,
                 coupling : float,
                 degeneracy : float,
                 theory : float):
        self.int2DScheme : str = None
        """ Scheme used to solve two-dimensional integrals
        allowed options include:
        
          - full: the inner integral is evaluated at arbitrary points 
	    selected automatically by the quadrature rule
        
	  - segregated: the inner integral is evaluated on a fixed 
	    grid that depends on the integrand that is being processed
        
        Segregated is usually faster than full but it could become 
	less accurate if the fixed points are not chosen correctly
        """
        self.threads : int = None
        """ Number of OMP threads for parallel calculations"""

    def print() -> None:
        """ Prints the content of the input structure"""
        pass

    def isEqual(inputs : qupled.Input) -> bool:
        """ Compares two input structures

        Returns:
            true if the two input structures are identical
        """
        pass

    
class SlfcGuess():
    """Class used to define an initial guess for the classical schemes (STLS, STLS-IET)"""
    def __init__(self):
        self.wvg : np.ndarray = None
        """ The wave-vector grid """
        self.slfc : np.ndarray = None
        """ The static local field correction """

        
class StlsInput(Input):
    """Class to handle the inputs related to the classical schemes.

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.
        theory: theory to solve
    """
    def __init__(self,
                 coupling : float,
                 degeneracy : float,
                 theory : float):
        self.chemicalPotential : list[float] = None
        """ Initial guess for the chemical potential """
        self.error : float = None
        """ minimum error for convergence """
        self.mixing : float = None
        """ mixing paramter """
        self.iet : str = None
        """ Classical-to-quantum mapping used in the iet schemes
        allowed options include:
        
          - standard: inversely proportional to the degeneracy parameter
        
	  - sqrt: inversely proportional to the square root of the sum
            of the squares of one and the degeneracy parameter

          - linear: inversely proportional to one plus the degeneracy
            parameter.
        
        Far from the ground state all mappings lead identical results, but at
        the ground state they can differ significantly (the standard
        mapping diverges)
        """
        self.matsubara : int = None
        """ NUmber of matsubara frequencies"""
        self.iterations : int = None
        """ Maximum number of iterations """
        self.outputFrequency : int = None
        """ Output frequency to write the recovery file """
        self.recoveryFile : str = None
        """ Name of the recovery file used to restart an unfinished simulation """
        self.guess : qupled.SlfcGuess = None
        """ Initial guess """
        self.resolution : float = None
        """ Resolution of the wave-vector grid """
        self.cutoff : float = None
        """ cutoff for the wave-vector grid """

class FreeEnergyIntegrand():
    """Class used to store the precomputed values of the free energy integrand"""
    def __init__(self):
        self.grid : np.ndarray = None
        """ The coupling parameter grid used to compute the free energy """
        self.integrand : np.ndarray = None
        """ The free energy integrand for various coupling parameter values """

        
class VSStlsInput(StlsInput):
    """Class to handle the inputs related to the classical VS-STLS scheme.

    Args:
        coupling: Coupling parameter.
        degeneracy: Degeneracy parameter.
    """
    def __init__(self,
                 coupling : float,
                 degeneracy : float):
        self.theory : str = None
        """ Name of the theory that is solved """
        self.alpha : float = None
        """ Initial guess for the free parameter """
        self.couplingResolution : float = None
        """ Resolution of the coupling parameter grid """
        self.degeneracyResolution : float = None
        """ Resolution of the degeneracy parameter grid """
        self.mixingAlpha : float = None
        """ Mixing parameter for the iterations to determine the free parameter """
        self.errorAlpha : float = None
        """ minimum error for convergence in the free parameter """
        self.freeEnergyIntegrand : qupled.FreeEnergyIntegrand = None
        """ Pre-computed free energy integrand """
    
class QstlsGuess():
    """Class used to define an initial guess for the quantum schemes (QSTLS, QSTLS-IET)"""
    def __init__(self):
        self.wvg : np.ndarray = None
        """ The wave-vector grid """
        self.ssf : np.ndarray = None
        """ The static structure factor """
        self.adr : np.ndarray = None
        """ The auxiliary density response """
        self.matsubara : int = None
        """ The number of matsubara frequencies """

        
class QstlsInput():
    """Class to handle the inputs related to the quantum schemes. """
    def __init__(self):
        self.guess : qupled.QstlsGuess = None
        """ Initial guess """
        self.fixed : float = None
        """ name of the file storing the fixed component of the auxiliary density 
	response in the QSTLS schmeme. Note: The QSTLS auxiliary density response 
	is computed also when the QSTLS-IET schemes are solved. So, if possible, it 
	is a good idea to set this property also when solving the QSTLS-IET scheme 
	because it allows to save quite some computational tim """
        self.fixediet : float = None
        """ name of the zip file storing the fixed components of the auxiliary density
	response in the QSTLS-IET schemes. Note: Whenever possible, it
	is a good idea to set this property when solving the QSTLS-IET schemes """

    def print() -> None:
        """ Prints the content of the input structure"""
        pass

    def isEqual(inputs : qupled.Input) -> bool:
        """ Compares two input structures

        Returns:
            true if the two input structures are identical
        """
        pass

    
class Stls():
    """Class to solve the classical schemes (STLS, STLS-IET)  

    Args:
        inputs: Input parameters.
    """
    def __init__(self,
                 inputs : qupled.StlsInput):
        self.bf : np.ndarray = None
        """ The bridge function adder (applies only to the IET schemes). Read only. """
        self.idr : np.ndarray = None
        """ The ideal density response. Read only. """
        self.recovery : str = None
        """ The name of the recovery file. Read only """
        self.sdr : np.ndarray = None
        """ The static density response. Read only. """
        self.slfc : np.ndarray = None
        """ The static local field correction. Read only. """
        self.ssf : np.ndarray = None
        """ The static structure factor. Read only. """
        self.ssfHF : np.ndarray = None
        """ The Hartree-Fock static structure factor. Read only. """
        self.wvg : np.ndarray = None
        """ The wave-vector grid """

    def compute() -> None:
        """ Solves the scheme """
        pass

    def getRdf() -> np.ndarray:
        """ Computes the radial distribution function

        Returns:
            The radial distribution function
        """
        pass


class VSStls(Stls):
    """Class to solve the classical VS-STLS scheme  

    Args:
        inputs: Input parameters.
    """
    def __init__(self,
                 inputs : qupled.VSStlsInput):
        self.freeEnergyIntegrand : np.ndarray = None
        """ Free energy integrand """
        self.freeEnergyGrid : np.ndarray = None
        """ Coupling parameter grid used to integrate the free energy integrand """
    
class Qstls(Stls):
    """Class to solve the quantum schemes (QSTLS, QSTLS-IET). This class inherits all the methods of its parent class.
    
    Args:
        inputs: Input parameters.
        qinputs: Input parameters (specific for quantum schemes)
    """
    def __init__(self,
                 inputs : qupled.StlsInput,
                 qinputs: qupled.QstlsInput):
        self.adr : np.ndarray = None
        """ The auxiliary density response. Read only. """

    
def computeRdf(rdfGrid : np.ndarray,
               wvg : np.ndarray,
               ssf : np.ndarray) -> np.ndarray:
    """ Computes the radial distribution function from the static structure factor.

    Args:
        rdfGrid: Grid used to compute the radial distribution function.
        wvg: The wave-vector grid.
        ssf: The static structure factor.
    
    Returns:
            The radial distribution function
    """
    pass

def computeInternalEnergy(wvg : np.ndarray,
                          ssf : np.ndarray,
                          coupling : float) -> float:
    """ Computes the internal energy from the static structure factor.

    Args:
        wvg: The wave-vector grid.
        ssf: The static structure factor.
        coupling: The coupling parameter.
    
    Returns:
        The internal energy
    """
    pass

def computeFreeEnergy(rsGrid: np.ndarray,
                      fxci : np.ndarray,
                      coupling : float) -> float:
    """ Computes the free energy for a given coupling parameter.

    Args:
        rsGrid: The coupling parameter grid.
        fxci: The free energy integrand.
        coupling: The coupling parameter.
    
    Returns:
        The free energy
    """
    pass
