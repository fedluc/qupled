from __future__ import annotations

from . import hf
from . import native


class ESA(hf.HF):
    """
    Class used to solve the ESA scheme.
    """

    def __init__(self):
        super().__init__()
        self.results: hf.Result = hf.Result()
        # Undocumented properties
        self.native_scheme_cls = native.ESA


class Input(hf.Input):
    """
    Class used to manage the input for the :obj:`qupled.esa.ESA` class.
    """

    def __init__(self, coupling: float, degeneracy: float):
        super().__init__(coupling, degeneracy)
        # Undocumented default values
        self.theory = "ESA"
