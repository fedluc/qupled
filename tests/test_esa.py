import qupled.esa as esa
import qupled.hf as hf
import qupled.native as native


def test_esa_inheritance():
    assert issubclass(esa.Solver, hf.Solver)


def test_esa_initialization(mocker):
    super_init = mocker.patch("qupled.hf.Solver.__init__")
    scheme = esa.Solver()
    super_init.assert_called_once()
    assert scheme.native_scheme_cls == native.ESA


def test_esa_input_inheritance():
    assert issubclass(esa.Input, hf.Input)


def test_esa_input_initialization(mocker):
    input = esa.Input(mocker.ANY, mocker.ANY)
    assert input.theory == "ESA"
