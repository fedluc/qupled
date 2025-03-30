import qupled.native as native
import qupled.esa as esa
import qupled.rpa as rpa


def test_esa_inheritance():
    assert issubclass(esa.ESA, rpa.Rpa)


def test_esa_initialization(mocker):
    super_init = mocker.patch("qupled.rpa.Rpa.__init__")
    scheme = esa.ESA()
    super_init.assert_called_once()
    assert scheme.native_scheme_cls == native.ESA


def test_esa_input_inheritance():
    assert issubclass(esa.Input, rpa.Input)


def test_esa_input_initialization(mocker):
    super_init = mocker.patch("qupled.rpa.Input.__init__")
    coupling = 1.5
    degeneracy = 3.0
    input = esa.Input(coupling, degeneracy)
    super_init.assert_called_once_with(coupling, degeneracy)
    assert input.theory == "ESA"
