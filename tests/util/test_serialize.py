from dataclasses import field

import numpy as np
import pytest

from qupled.util import serialize


@serialize.serializable_dataclass
class Inner:
    value: int


@serialize.serializable_dataclass
class Outer:
    required: int
    inner: Inner = field(default_factory=lambda: Inner(1))
    array: np.ndarray = None
    values: list[int] = field(default_factory=list)


@serialize.serializable_dataclass
class Validated:
    value: int

    def __post_init__(self):
        if self.value < 0:
            raise ValueError("value must be non-negative")


@pytest.mark.unit
def test_from_dict_preserves_dataclass_defaults():
    result = Outer.from_dict({"required": 10})

    assert result.required == 10
    assert result.inner == Inner(1)
    assert result.array is None
    assert result.values == []


@pytest.mark.unit
def test_from_dict_recursively_converts_values_before_init():
    result = Outer.from_dict(
        {
            "required": 10,
            "inner": {"value": 2},
            "array": [1.0, 2.0],
            "values": [3, 4],
        }
    )

    assert result.inner == Inner(2)
    np.testing.assert_array_equal(result.array, np.array([1.0, 2.0]))
    assert result.values == [3, 4]


@pytest.mark.unit
def test_from_dict_uses_dataclass_required_field_validation():
    with pytest.raises(TypeError):
        Outer.from_dict({})


@pytest.mark.unit
def test_from_dict_uses_post_init_validation():
    with pytest.raises(ValueError, match="value must be non-negative"):
        Validated.from_dict({"value": -1})
