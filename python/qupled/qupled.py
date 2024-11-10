# Placeholders used to document the c++ classes exposed with python::boost
# NOTE: This is just a workaround, a more stable solution will be developed in the future

from __future__ import annotations


class Rpa:
    pass


class ESA:
    pass


class Stls:
    pass


class VSStls:
    pass


class Qstls:
    pass


class QVSStls:
    pass


class RpaInput:
    pass


class IterationInput:
    pass


class StlsInput(RpaInput, IterationInput):
    pass


class VSInput:
    pass


class VSStlsInput(StlsInput, VSInput):
    pass


class QstlsInput(RpaInput, IterationInput):
    pass


class QVSStlsInput(QstlsInput, VSInput):
    pass


class StlsGuess:
    """Class used to define an initial guess for the classical schemes (STLS, STLS-IET)."""

    def __init__(self):
        self.wvg: np.ndarray = None
        """ The wave-vector grid """
        self.slfc: np.ndarray = None
        """ The static local field correction """


class QstlsGuess:
    """Class used to define an initial guess for the quantum schemes (QSTLS, QSTLS-IET)."""

    def __init__(self):
        self.wvg: np.ndarray = None
        """ The wave-vector grid """
        self.ssf: np.ndarray = None
        """ The static structure factor """
        self.adr: np.ndarray = None
        """ The auxiliary density response """
        self.matsubara: int = None
        """ The number of matsubara frequencies """


class FreeEnergyIntegrand:
    """Class used to store the precomputed values of the free energy integrand for the VS schemes."""

    def __init__(self):
        self.grid: np.ndarray = None
        """ The coupling parameter grid used to compute the free energy """
        self.integrand: np.ndarray = None
        """ The free energy integrand for various coupling parameter values """
        self.alpha: np.ndarray = None
