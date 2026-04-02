import pytest

from qupled import native
from qupled.schemes import esa, hf, rpa


@pytest.mark.unit
def test_esa_inheritance():
    assert issubclass(esa.Solver, hf.Solver)


@pytest.mark.unit
def test_esa_initialization(mocker):
    super_init = mocker.patch("qupled.schemes.hf.Solver.__init__")
    scheme = esa.Solver()
    super_init.assert_called_once()
    assert scheme.native_scheme_cls == native.ESA


@pytest.mark.unit
def test_esa_input_inheritance():
    assert issubclass(esa.Input, hf.Input)


@pytest.mark.unit
def test_esa_input_initialization(mocker):
    input = esa.Input(mocker.ANY, mocker.ANY)
    assert input.theory == "ESA"
