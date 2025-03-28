from qupled.base import QuantumIterativeScheme
from qupled.native import Qstls as NativeQstls
from qupled.qstls import Qstls
from qupled.stls import Stls


def test_qstls_inheritance():
    assert issubclass(Qstls, QuantumIterativeScheme)


def test_qstls_compute(mocker):
    super_compute = mocker.patch("qupled.base.QuantumIterativeScheme.compute")
    QstlsInput = mocker.patch("qupled.native.QstlsInput")
    qstls_input = Qstls.Input(coupling=1.0, degeneracy=2.0)
    native_input = mocker.MagicMock()
    QstlsInput.return_value = native_input
    qstls = Qstls()
    result = mocker.MagicMock()
    mocker.patch.object(qstls, "Result", return_value=result)
    qstls.compute(qstls_input)
    super_compute.assert_called_once_with(
        qstls_input, NativeQstls, native_input, result
    )


def test_qstls_input_inheritance():
    assert issubclass(Qstls.Input, Stls.Input)


def test_qstls_input_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    qstls_input = Qstls.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert qstls_input.theory == "QSTLS"


def test_qstls_result_inheritance():
    assert issubclass(Qstls.Result, Stls.Result)


def test_qstls_result_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.Result.__init__")
    qstls_result = Qstls.Result()
    super_init.assert_called_once()
    assert qstls_result.adr is None
