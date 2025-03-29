import pytest
from qupled.native import Rpa as NativeRpa
from qupled.base import ClassicScheme, Input
from qupled.rpa import Rpa


def test_esa_inheritance():
    assert issubclass(Rpa, ClassicScheme)


def test_rpa_compute(mocker):
    super_compute = mocker.patch("qupled.base.ClassicScheme.compute")
    input = mocker.ANY
    scheme = Rpa()
    scheme.compute(input)
    super_compute.assert_called_once_with(input, NativeRpa, mocker.ANY, mocker.ANY)


def test_esa_input_inheritance():
    assert issubclass(Rpa.Input, Input)


def test_rpa_input_initialization(mocker):
    super_init = mocker.patch("qupled.base.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    rpa_input = Rpa.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert rpa_input.theory == "RPA"
