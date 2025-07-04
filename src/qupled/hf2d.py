from __future__ import annotations

from . import hf
from . import native
from . import serialize


class Solver(hf.Solver):
    """
    Class used to solve the HF scheme in 2D or 3D.
    """
    native_scheme_cls = native.HF  # Using same class for both dimensions

@serialize.serializable_dataclass
class Input(hf.Input):
    """
    Class used to manage the input for HF calculations.
    """
    dimension: str = "2D" 
    theory: str = "HF"     

    def __post_init__(self):
        super().__post_init__()
        self.dimension = "2D"  
