import pytest

import qupled.qstls as qstls
import qupled.qstlsiet as qstlsiet
import qupled.stlsiet as stlsiet


@pytest.fixture
def input(mocker):
    return mocker.Mock()


@pytest.fixture
def scheme():
    return qstlsiet.Solver()


def test_qstls_iet_inheritance():
    assert issubclass(qstlsiet.Solver, qstls.Solver)


def test_qstls_iet_initialization(mocker):
    super_init = mocker.patch("qupled.qstls.Solver.__init__")
    scheme = qstlsiet.Solver()
    super_init.assert_called_once()
    assert isinstance(scheme.results, stlsiet.Result)


def test_qstls_iet_input_inheritance():
    assert issubclass(qstlsiet.Input, (stlsiet.Input, qstls.Input))


def test_qstls_iet_input_initialization_valid_theory(mocker):
    theory = "QSTLS-HNC"
    input = qstlsiet.Input(mocker.ANY, mocker.ANY, theory=theory)
    assert input.theory == theory


def test_qstls_iet_input_initialization_invalid_theory():
    with pytest.raises(ValueError):
        qstlsiet.Input(1.0, 1.0, "INVALID-THEORY")


def test_qstls_iet_input_initialization_default_theory():
    with pytest.raises(ValueError):
        qstlsiet.Input(1.0, 1.0)
