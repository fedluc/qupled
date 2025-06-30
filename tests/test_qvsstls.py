import pytest

import qupled.qstls as qstls
import qupled.qvsstls as qvsstls
import qupled.vsstls as vsstls


@pytest.fixture
def input(mocker):
    return mocker.Mock()


@pytest.fixture
def scheme():
    return qvsstls.Solver()


def test_qvsstls_inheritance():
    assert issubclass(qvsstls.Solver, vsstls.Solver)


def test_compute(mocker):
    find_fixed_adr_in_database = mocker.patch(
        "qupled.qstls.Solver.find_fixed_adr_in_database"
    )
    super_compute = mocker.patch("qupled.vsstls.Solver.compute")
    inputs = mocker.ANY
    scheme = qvsstls.Solver()
    scheme.compute(inputs)
    find_fixed_adr_in_database.assert_called_once_with(scheme, inputs)
    super_compute.assert_called_once_with(inputs)


def test_get_free_energy_integrand(mocker):
    run_id = mocker.ANY
    database_name = mocker.ANY
    get_free_energy_integrand = mocker.patch(
        "qupled.vsstls.Solver.get_free_energy_integrand"
    )
    result = qvsstls.Solver.get_free_energy_integrand(run_id, database_name)
    get_free_energy_integrand.assert_called_once_with(run_id, database_name)
    assert result == get_free_energy_integrand.return_value


def test_qvsstls_input_inheritance():
    assert issubclass(qvsstls.Input, (qstls.Input, vsstls.Input))


def test_qvsstls_input_initialization_valid_theory(mocker):
    input = qvsstls.Input(mocker.ANY, mocker.ANY)
    assert input.theory == "QVSSTLS"
