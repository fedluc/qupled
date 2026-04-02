import pytest
import numpy as np

from qupled import native
from qupled.schemes import hf, rpa


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
    super_init = mocker.patch("qupled.schemes.hf.Solver.__init__")
    scheme = rpa.Solver()
    super_init.assert_called_once()
    assert scheme.native_scheme_cls == native.Rpa


def test_rpa_input_inheritance():
    assert issubclass(rpa.Input, hf.Input)


def test_rpa_input_initialization(mocker):
    input = rpa.Input(mocker.ANY, mocker.ANY)
    assert input.theory == "RPA"


def test_rpa_result_inheritance():
    assert issubclass(rpa.Result, hf.Result)


def test_rpa_result_compute_itcf(mocker, results, inputs):
    native_compute_itcf = mocker.patch("qupled.native.compute_itcf")
    native_input = mocker.Mock()
    native_inputs_cls = mocker.patch.object(native, "Input", return_value=native_input)
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.lfc = np.array([4.0, 5.0, 6.0])
    results.chemical_potential = 0.5
    results.idr = np.array([7.0, 8.0, 9.0])
    custom_tau = np.array([0.0, 0.2, 0.4])
    native_compute_itcf.return_value = np.array([[13.0, 14.0, 15.0]])
    results.compute_itcf(inputs, custom_tau)
    assert np.allclose(results.tau, custom_tau)
    native_compute_itcf.assert_called_once_with(
        native_input,
        results.wvg,
        custom_tau,
        results.chemical_potential,
        results.idr,
        results.lfc,
    )
    assert np.allclose(results.itcf, np.array([[13.0, 14.0, 15.0]]))
