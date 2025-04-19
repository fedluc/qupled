from __future__ import annotations

from . import qstls
from . import stlsiet


class QstlsIet(qstls.Qstls):
    """
    Class used to solve the Qstls-IET schemes.
    """

    def __init__(self):
        super().__init__()
        self.results: Result = Result()


# Input class
class Input(stlsiet.Input, qstls.Input):
    """
    Class used to manage the input for the :obj:`qupled.qstlsiet.QStlsIet` class.
    Accepted theories: ``QSTLS-HNC``, ``QSTLS-IOI`` and ``QSTLS-LCT``.
    """

    def __init__(self, coupling: float, degeneracy: float, theory: str):
        stlsiet.Input.__init__(self, coupling, degeneracy, "STLS-HNC")
        qstls.Input.__init__(self, coupling, degeneracy)
        if theory not in {"QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"}:
            raise ValueError("Invalid dielectric theory")
        self.theory = theory
        self.fixed_iet = ""
        """
        Name of the zip file storing the iet part of the fixed components
        of the auxiliary density response. Default = ``""``
        """


# Result class
class Result(stlsiet.Result, qstls.Result):
    """
    Class used to store the results for the :obj:`qupled.qstlsiet.QstlsIet` class.
    """

    def __init__(self):
        super().__init__()
