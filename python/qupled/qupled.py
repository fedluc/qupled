# Placeholders used to document the c++ classes exposed with python::boost

from __future__ import annotations

class RpaInput():
    """Class to handle the inputs related to the classical RPA scheme
    and the hybrid ESA scheme.
    """
    def __init__(self):
        self.coupling : float
        """ Coupling parameter """
        self.degeneracy : float
        """ Degeneracy parameter """
        self.theory : float
        """ Theory to solve """
        self.intError : float = None
        """ Accuracy (expressed as a relative error) in the computation of the integrals """
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
        self.chemicalPotential : list[float] = None
        """ Initial guess for the chemical potential """
        self.matsubara : int = None
        """ Number of matsubara frequencies"""
        self.resolution : float = None
        """ Resolution of the wave-vector grid """
        self.cutoff : float = None
        """ cutoff for the wave-vector grid """

class StlsInput(RpaInput):
    """Class to handle the inputs related to the classical STLS and STLS-IET schemes."""
    def __init__(self):
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
        self.iterations : int = None
        """ Maximum number of iterations """
        self.outputFrequency : int = None
        """ Output frequency to write the recovery file """
        self.recoveryFile : str = None
        """ Initial guess """
        
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
        self.alpha : list[float] = None
        """ Initial guess for the free parameter """
        self.couplingResolution : float = None
        """ Resolution of the coupling parameter grid """
        self.degeneracyResolution : float = None
        """ Resolution of the degeneracy parameter grid """
        self.errorAlpha : float = None
        """ Minimum error for convergence in the free parameter """
        self.iterationsAlpha : int = None
        """ Maximum number of iterations to determine the free parameter """
        self.freeEnergyIntegrand : qupled.FreeEnergyIntegrand = None
        """ Pre-computed free energy integrand """

        
class QstlsInput(StlsInput):
    """Class to handle the inputs related to the quantum schemes. """
    def __init__(self):
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
        
class QVSStlsInput(VSStlsInput, QstlsInput):
    """Class to handle the inputs related to the quantum VS-STLS scheme."""
    pass


class SlfcGuess():
    """Class used to define an initial guess for the classical schemes (STLS, STLS-IET)."""
    def __init__(self):
        self.wvg : np.ndarray = None
        """ The wave-vector grid """
        self.slfc : np.ndarray = None
        """ The static local field correction """

class QstlsGuess():
    """Class used to define an initial guess for the quantum schemes (QSTLS, QSTLS-IET)."""
    def __init__(self):
        self.wvg : np.ndarray = None
        """ The wave-vector grid """
        self.ssf : np.ndarray = None
        """ The static structure factor """
        self.adr : np.ndarray = None
        """ The auxiliary density response """
        self.matsubara : int = None
        """ The number of matsubara frequencies """
        
class FreeEnergyIntegrand():
    """Class used to store the precomputed values of the free energy integrand for the VS-STLS scheme."""
    def __init__(self):
        self.grid : np.ndarray = None
        """ The coupling parameter grid used to compute the free energy """
        self.integrand : np.ndarray = None
        """ The free energy integrand for various coupling parameter values """
        self.alpha : np.ndarray = None
