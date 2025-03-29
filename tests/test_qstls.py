from qupled.base import QuantumIterativeScheme
from qupled.native import Qstls as NativeQstls
from qupled.qstls import Qstls
from qupled.stls import Stls


def test_qstls_inheritance():
    assert issubclass(Qstls, QuantumIterativeScheme)


def test_qstls_compute(mocker):
    super_compute = mocker.patch("qupled.base.ClassicScheme.compute")
    input = mocker.ANY
    scheme = Qstls()
    scheme.compute(input)
    super_compute.assert_called_once_with(input, NativeQstls, mocker.ANY, mocker.ANY)


def test_qstls_input_inheritance():
    assert issubclass(Qstls.Input, Stls.Input)


def test_qstls_input_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    input = Qstls.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "QSTLS"


def test_qstls_result_inheritance():
    assert issubclass(Qstls.Result, Stls.Result)


def test_qstls_result_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.Result.__init__")
    result = Qstls.Result()
    super_init.assert_called_once()
    assert result.adr is None
