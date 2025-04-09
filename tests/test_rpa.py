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
    scheme = rpa.Rpa()
    scheme.db_handler = mocker.Mock()
    return scheme


def test_rpa_initialization(mocker):
    super_init = mocker.patch("qupled.hf.HF.__init__")
    scheme = rpa.Rpa()
    super_init.assert_called_once()
    assert scheme.native_scheme_cls == native.Rpa


def test_rpa_input_inheritance():
    assert issubclass(rpa.Input, hf.Input)


def test_rpa_input_initialization(mocker):
    super_init = mocker.patch("qupled.hf.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    input = rpa.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "RPA"
