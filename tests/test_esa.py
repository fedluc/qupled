import qupled.esa as esa
import qupled.hf as hf
import qupled.native as native


def test_esa_inheritance():
    assert issubclass(esa.ESA, hf.HF)


def test_esa_initialization(mocker):
    super_init = mocker.patch("qupled.hf.HF.__init__")
    scheme = esa.ESA()
    super_init.assert_called_once()
    assert scheme.native_scheme_cls == native.ESA


def test_esa_input_inheritance():
    assert issubclass(esa.Input, hf.Input)


def test_esa_input_initialization(mocker):
    super_init = mocker.patch("qupled.hf.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    input = esa.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "ESA"
