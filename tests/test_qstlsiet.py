import pytest

import qupled.qstls as qstls
import qupled.qstlsiet as qstlsiet
import qupled.stlsiet as stlsiet


@pytest.fixture
def input(mocker):
    return mocker.Mock()


@pytest.fixture
def scheme():
    return qstlsiet.QstlsIet()


def test_qstls_iet_inheritance():
    assert issubclass(qstlsiet.QstlsIet, qstls.Qstls)


def test_qstls_iet_initialization(mocker):
    super_init = mocker.patch("qupled.qstls.Qstls.__init__")
    scheme = qstlsiet.QstlsIet()
    super_init.assert_called_once()
    assert isinstance(scheme.results, stlsiet.Result)


def test_qstls_iet_input_inheritance():
    assert issubclass(qstlsiet.Input, (stlsiet.Input, qstls.Input))


def test_qstls_iet_input_initialization_valid_theory(mocker):
    qstls_init = mocker.patch("qupled.qstls.Input.__init__")
    stls_iet_init = mocker.patch("qupled.stlsiet.Input.__init__")
    coupling = 1.0
    degeneracy = 1.0
    theory = "QSTLS-HNC"
    input = qstlsiet.Input(coupling, degeneracy, theory)
    # qstls_init.assert_called_once_with(input, coupling, degeneracy)
    # stls_iet_init.assert_called_once_with(input, coupling, degeneracy, "STLS-HNC")
    # assert input.theory == theory


def test_qstls_iet_input_initialization_invalid_theory():
    with pytest.raises(ValueError):
        qstlsiet.Input(1.0, 1.0, "INVALID-THEORY")
