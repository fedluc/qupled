import pytest
from unittest import mock
from qupled.esa import ESA
from qupled.native import ESA as NativeESA
from qupled.base import ClassicScheme
from qupled.rpa import Rpa

def test_esa_inheritance():
    assert issubclass(ESA, ClassicScheme)

@mock.patch("qupled.base.ClassicScheme.compute")
@mock.patch("qupled.base.Result")
@mock.patch("qupled.native.RpaInput")
def test_esa_compute(RpaInput, Result, super_compute, esa_input):
    esa_input = ESA.Input(coupling=1.0, degeneracy=2.0)
    native_input = mock.MagicMock()
    result = mock.MagicMock()
    RpaInput.return_value = native_input
    Result.return_value = result
    esa = ESA()
    esa.compute(esa_input)
    super_compute.assert_called_once_with(esa_input, NativeESA, native_input, result)


def test_esa_input_inheritance():
    assert issubclass(ESA.Input, Rpa.Input)

@mock.patch("qupled.rpa.Rpa.Input.__init__")
def test_esa_input_initialization(super_init):
    coupling = 1.5
    degeneracy = 3.0
    esa_input = ESA.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert esa_input.theory == "ESA"

