from __future__ import annotations

from . import native
from . import qstls
from . import stlsiet


class QstlsIet(qstls.Qstls):
    """
    Class used to solve the Qstls-IET schemes.
    """

    def __init__(self):
        super().__init__()
        self.results: stlsiet.Result = stlsiet.Result()
        self.native_scheme_cls = native.QstlsIet
        self.native_inputs_cls = native.QstlsIetInput


# Input class
class Input(stlsiet.Input, qstls.Input):
    """
    Class used to manage the input for the :obj:`qupled.qstlsiet.QStlsIet` class.
    Accepted theories: ``QSTLS-HNC``, ``QSTLS-IOI`` and ``QSTLS-LCT``.
    """

    def __init__(self, coupling: float, degeneracy: float, theory: str):
        super().__init__(coupling, degeneracy, "STLS-HNC")
        if theory not in {"QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"}:
            raise ValueError("Invalid dielectric theory")
        self.theory = theory
