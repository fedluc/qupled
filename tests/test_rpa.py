import pytest

import qupled.hf as hf
import qupled.native as native
import qupled.rpa as rpa


@pytest.fixture
def inputs():
    return rpa.Input(coupling=1.0, degeneracy=2.0)


@pytest.fixture
def results():
    return rpa.Result()


@pytest.fixture
def scheme(mocker):
    scheme = rpa.Solver()
    scheme.db_handler = mocker.Mock()
    return scheme


def test_rpa_initialization(mocker):
    super_init = mocker.patch("qupled.hf.Solver.__init__")
    scheme = rpa.Solver()
    super_init.assert_called_once()
    assert scheme.native_scheme_cls == native.Rpa


def test_rpa_input_inheritance():
    assert issubclass(rpa.Input, hf.Input)


def test_rpa_input_initialization(mocker):
    input = rpa.Input(mocker.ANY, mocker.ANY)
    assert input.theory == "RPA"
