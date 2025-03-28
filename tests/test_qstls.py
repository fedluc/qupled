from unittest import mock
from qupled.base import QuantumIterativeScheme
from qupled.native import Qstls as NativeQstls
from qupled.qstls import Qstls
from qupled.stls import Stls


def test_qstls_inheritance():
    assert issubclass(Qstls, QuantumIterativeScheme)


@mock.patch("qupled.base.ClassicScheme.compute")
@mock.patch("qupled.native.QstlsInput")
def test_qstls_compute(QstlsInput, super_compute, mocker):
    qstls_input = Qstls.Input(coupling=1.0, degeneracy=2.0)
    native_input = mock.MagicMock()
    QstlsInput.return_value = native_input
    qstls = Qstls()
    result = mock.MagicMock()
    mocker.patch.object(qstls, "Result", return_value=result)
    qstls.compute(qstls_input)
    super_compute.assert_called_once_with(
        qstls_input, NativeQstls, native_input, result
    )


def test_qstls_input_inheritance():
    assert issubclass(Qstls.Input, Stls.Input)


@mock.patch("qupled.stls.Stls.Input.__init__")
def test_qstls_input_initialization(super_init):
    coupling = 1.5
    degeneracy = 3.0
    qstls_input = Qstls.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert qstls_input.theory == "QSTLS"


def test_qstls_result_inheritance():
    assert issubclass(Qstls.Result, Stls.Result)


@mock.patch("qupled.stls.Stls.Result.__init__")
def test_qstls_result_initialization(super_init):
    qstls_input = Qstls.Result()
    super_init.assert_called_once()
    assert qstls_input.adr is None
