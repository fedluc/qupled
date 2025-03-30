from __future__ import annotations
from . import native
from . import rpa


class ESA(rpa.Rpa):

    def __init__(self):
        super().__init__()
        # Undocumented properties
        self.native_scheme_cls = native.ESA


class Input(rpa.Input):
    """
    Class used to manage the input for the :obj:`qupled.classic.ESA` class.
    """

    def __init__(self, coupling: float, degeneracy: float):
        super().__init__(coupling, degeneracy)
        # Undocumented default values
        self.theory = "ESA"
