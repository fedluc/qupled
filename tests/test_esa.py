from qupled.esa import ESA
from qupled.native import ESA as NativeESA
from qupled.base import ClassicScheme
from qupled.rpa import Rpa


def test_esa_inheritance():
    assert issubclass(ESA, ClassicScheme)


def test_esa_compute(mocker):
    super_compute = mocker.patch("qupled.base.ClassicScheme.compute")
    Result = mocker.patch("qupled.base.Result")
    RpaInput = mocker.patch("qupled.native.RpaInput")
    esa_input = ESA.Input(coupling=1.0, degeneracy=2.0)
    native_input = mocker.ANY
    result = mocker.ANY
    RpaInput.return_value = native_input
    Result.return_value = result
    esa = ESA()
    esa.compute(esa_input)
    super_compute.assert_called_once_with(esa_input, NativeESA, native_input, result)


def test_esa_input_inheritance():
    assert issubclass(ESA.Input, Rpa.Input)


def test_esa_input_initialization(mocker):
    super_init = mocker.patch("qupled.rpa.Rpa.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    esa_input = ESA.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert esa_input.theory == "ESA"
