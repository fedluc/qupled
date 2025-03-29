from qupled.rpa import Rpa
from qupled.stls import Stls
from qupled.native import Stls as NativeStls
from qupled.base import IterativeScheme, Result


def test_stls_inheritance():
    assert issubclass(Stls, IterativeScheme)


def test_stls_compute(mocker):
    super_compute = mocker.patch("qupled.base.ClassicScheme.compute")
    input = mocker.ANY
    scheme = Stls()
    scheme.compute(input)
    super_compute.assert_called_once_with(input, NativeStls, mocker.ANY, mocker.ANY)


def test_stls_input_inheritance():
    assert issubclass(Stls.Input, Rpa.Input)


def test_stls_input_initialization(mocker):
    super_init = mocker.patch("qupled.base.Input.__init__")
    guess = mocker.patch("qupled.stls.Stls.Guess")
    coupling = 1.5
    degeneracy = 3.0
    input = Stls.Input(coupling, degeneracy)
    assert input.error == 1.0e-5
    assert input.mixing == 1.0
    assert input.iterations == 1000
    assert input.output_frequency == 10
    assert input.recovery_file == ""
    assert input.guess == guess.return_value
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "STLS"


def test_stls_result_inheritance():
    assert issubclass(Stls.Result, Rpa.Result)


def test_stls_result_initialization(mocker):
    super_init = mocker.patch("qupled.base.Result.__init__")
    results = Stls.Result()
    assert results.error is None
    super_init.assert_called_once()
