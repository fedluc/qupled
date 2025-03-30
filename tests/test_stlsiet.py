import pytest
import qupled.stlsiet as stlsiet
import qupled.stls as stls


def test_stls_iet_inheritance():
    assert issubclass(stlsiet.StlsIet, stls.Stls)


def test_stls_iet_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Stls.__init__")
    scheme = stlsiet.StlsIet()
    super_init.assert_called_once()
    assert isinstance(scheme.results, stlsiet.Result)


def test_stls_iet_input_inheritance():
    assert issubclass(stlsiet.Input, stls.Input)


def test_stls_iet_input_initialization_valid_theory(mocker):
    super_init = mocker.patch("qupled.stls.Input.__init__")
    coupling = 1.0
    degeneracy = 1.0
    theory = "STLS-HNC"
    input = stlsiet.Input(coupling, degeneracy, theory)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == theory


def test_stls_iet_input_initialization_invalid_theory():
    with pytest.raises(ValueError):
        stlsiet.Input(1.0, 1.0, "INVALID-THEORY")


def test_stls_iet_result_inheritance():
    assert issubclass(stlsiet.Result, stls.Result)


def test_stls_iet_result_initialization(mocker):
    super_init = mocker.patch("qupled.stls.Result.__init__")
    results = stlsiet.Result()
    assert results.bf is None
    super_init.assert_called_once()
