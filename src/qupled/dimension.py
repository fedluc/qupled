from enum import Enum

from . import native


class Dimension(Enum):
    """
    An enumerator representing spatial dimensions.

    Values:

       ``_2D``: Represents a two-dimensional space.

       ``_3D``: Represents a three-dimensional space.
    """

    _2D = "D2"
    _3D = "D3"

    def to_dict(self):
        return {"dimension": self.value}

    @classmethod
    def from_dict(cls, d):
        value = d.get("dimension", None)
        if value is None:
            raise ValueError("No 'dimension' key found in dict")
        for dim in cls:
            if dim.value == value:
                return dim
        raise ValueError(f"Unknown dimension value: {value}")

    def to_native(self):
        return getattr(native.Dimension, self.value)
