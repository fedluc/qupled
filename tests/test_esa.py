from qupled.esa import ESA
from qupled.native import ESA as NativeESA
from qupled.base import ClassicScheme
from qupled.rpa import Rpa


def test_esa_inheritance():
    assert issubclass(ESA, ClassicScheme)


def test_esa_compute(mocker):
    super_compute = mocker.patch("qupled.base.ClassicScheme.compute")
    input = mocker.ANY
    scheme = ESA()
    scheme.compute(input)
    super_compute.assert_called_once_with(input, NativeESA, mocker.ANY, mocker.ANY)


def test_esa_input_inheritance():
    assert issubclass(ESA.Input, Rpa.Input)


def test_esa_input_initialization(mocker):
    super_init = mocker.patch("qupled.rpa.Rpa.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    input = ESA.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "ESA"
