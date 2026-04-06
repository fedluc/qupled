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


@pytest.mark.unit
def test_rpa_result_inheritance():
    assert issubclass(rpa.Result, hf.Result)


@pytest.mark.unit
def test_result_compute_itcf_with_default_grid(mocker, results, inputs):
    # inputs.degeneracy = 2.0, so default tau = np.arange(0.0, 1.1, 0.1) / 2.0
    native_compute_itcf = mocker.patch("qupled.native.compute_itcf")
    native_input = mocker.Mock()
    mocker.patch.object(native, "Input", return_value=native_input)
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.lfc = np.array([4.0, 5.0, 6.0])
    results.chemical_potential = 0.5
    results.idr = np.array([7.0, 8.0, 9.0])
    native_compute_itcf.return_value = np.array([[10.0, 11.0, 12.0]])
    results.compute_itcf(inputs)
    assert np.allclose(results.tau, np.arange(0.0, 1.1, 0.1) / inputs.degeneracy)
    native_compute_itcf.assert_called_once_with(
        native_input,
        results.wvg,
        results.tau,
        results.chemical_potential,
        results.idr,
        results.lfc,
    )
    assert np.allclose(results.itcf, np.array([[10.0, 11.0, 12.0]]))


@pytest.mark.unit
def test_rpa_result_compute_itcf_with_custom_grid(mocker, results, inputs):
    # inputs.degeneracy = 2.0; all custom_tau values are <= 1/2.0, so none filtered
    native_compute_itcf = mocker.patch("qupled.native.compute_itcf")
    native_input = mocker.Mock()
    mocker.patch.object(native, "Input", return_value=native_input)
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
        results.tau,
        results.chemical_potential,
        results.idr,
        results.lfc,
    )
    assert np.allclose(results.itcf, np.array([[13.0, 14.0, 15.0]]))


@pytest.mark.unit
def test_rpa_result_invoke_native_itcf(mocker, results, inputs):
    native_compute_itcf = mocker.patch("qupled.native.compute_itcf")
    native_input = mocker.Mock()
    mocker.patch.object(native, "Input", return_value=native_input)
    results.wvg = np.array([1.0, 2.0, 3.0])
    results.lfc = np.array([4.0, 5.0, 6.0])
    results.tau = np.array([0.0, 0.1, 0.2])
    results.chemical_potential = 0.5
    results.idr = np.array([7.0, 8.0, 9.0])
    expected = np.array([[10.0, 11.0, 12.0]])
    native_compute_itcf.return_value = expected
    result = results._invoke_native_itcf(inputs)
    native_compute_itcf.assert_called_once_with(
        native_input,
        results.wvg,
        results.tau,
        results.chemical_potential,
        results.idr,
        results.lfc,
    )
    assert np.allclose(result, expected)
