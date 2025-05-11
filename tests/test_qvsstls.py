import pytest

import qupled.qstls as qstls
import qupled.qvsstls as qvsstls
import qupled.vsstls as vsstls


@pytest.fixture
def input(mocker):
    return mocker.Mock()


@pytest.fixture
def scheme():
    return qvsstls.QVSStls()


def test_qvsstls_inheritance():
    assert issubclass(qvsstls.QVSStls, vsstls.VSStls)


def test_compute(mocker):
    find_fixed_adr_in_database = mocker.patch(
        "qupled.qstls.Qstls.find_fixed_adr_in_database"
    )
    super_compute = mocker.patch("qupled.vsstls.VSStls.compute")
    inputs = mocker.ANY
    scheme = qvsstls.QVSStls()
    scheme.compute(inputs)
    find_fixed_adr_in_database.assert_called_once_with(scheme, inputs)
    super_compute.assert_called_once_with(inputs)


def test_get_free_energy_integrand(mocker):
    run_id = mocker.ANY
    database_name = mocker.ANY
    get_free_energy_integrand = mocker.patch(
        "qupled.vsstls.VSStls.get_free_energy_integrand"
    )
    result = qvsstls.QVSStls.get_free_energy_integrand(run_id, database_name)
    get_free_energy_integrand.assert_called_once_with(run_id, database_name)
    assert result == get_free_energy_integrand.return_value


def test_qvsstls_input_inheritance():
    assert issubclass(qvsstls.Input, (qstls.Input, vsstls.Input))


def test_qvsstls_input_initialization_valid_theory(mocker):
    qstls_init = mocker.patch("qupled.qstls.Input.__init__")
    vsstls_init = mocker.patch("qupled.vsstls.Input.__init__")
    coupling = 1.0
    degeneracy = 1.0
    input = qvsstls.Input(coupling, degeneracy)
    qstls_init.assert_called_once_with(input, coupling, degeneracy)
    vsstls_init.assert_called_once_with(input, coupling, degeneracy)
    assert input.theory == "QVSSTLS"
