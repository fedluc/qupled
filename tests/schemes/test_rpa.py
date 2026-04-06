import pytest

from qupled import native
from qupled.schemes import hf, rpa


@pytest.fixture
def inputs():
    return rpa.Input(coupling=1.0, degeneracy=2.0)


@pytest.fixture
def scheme(mocker):
    scheme = rpa.Solver()
    scheme.db_handler = mocker.Mock()
    return scheme


@pytest.mark.unit
def test_rpa_initialization(mocker):
    super_init = mocker.patch("qupled.schemes.hf.Solver.__init__")
    scheme = rpa.Solver()
    super_init.assert_called_once()
    assert scheme.native_scheme_cls == native.Rpa


@pytest.mark.unit
def test_rpa_input_inheritance():
    assert issubclass(rpa.Input, hf.Input)


@pytest.mark.unit
def test_rpa_input_initialization(mocker):
    input = rpa.Input(mocker.ANY, mocker.ANY)
    assert input.theory == "RPA"
