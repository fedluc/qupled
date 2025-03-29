import pytest
from qupled.stlsiet import StlsIet
from qupled.stls import Stls
from qupled.native import Stls as NativeStls
from qupled.base import IterativeScheme, Input, Result


def test_stls_iet_inheritance():
    assert issubclass(StlsIet, IterativeScheme)


def test_stls_iet_compute(mocker):
    super_compute = mocker.patch("qupled.base.ClassicScheme.compute")
    input = mocker.ANY
    scheme = StlsIet()
    scheme.compute(input)
    super_compute.assert_called_once_with(input, NativeStls, mocker.ANY, mocker.ANY)


def test_stls_iet_input_inheritance():
    assert issubclass(StlsIet.Input, Stls.Input)


def test_stls_iet_input_initialization_valid_theory(mocker):
    super_init = mocker.patch("qupled.stls.Stls.Input.__init__")
    coupling = 1.0
    degeneracy = 1.0
    theory = "STLS-HNC"
    input = StlsIet.Input(coupling, degeneracy, theory)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == theory


def test_stls_iet_input_initialization_invalid_theory():
    with pytest.raises(ValueError):
        StlsIet.Input(1.0, 1.0, "INVALID-THEORY")


def test_stls_iet_result_inheritance():
    assert issubclass(StlsIet.Result, Stls.Result)


def test_stls_iet_result_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.Result.__init__")
    results = StlsIet.Result()
    assert results.bf is None
    super_init.assert_called_once()
