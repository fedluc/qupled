import importlib
import os
import shutil
import sys

import matplotlib.pyplot as plt
import pytest

from qupled.database.database_handler import DataBaseHandler, DATABASE_DIRECTORY


@pytest.fixture(autouse=True)
def run_before_each_test():
    examples_dir = os.path.abspath("docs")
    if examples_dir not in sys.path:
        sys.path.insert(0, examples_dir)
    yield
    if examples_dir in sys.path:
        sys.path.remove(examples_dir)


@pytest.fixture(autouse=True)
def run_after_each_test():
    yield
    output_dir = DATABASE_DIRECTORY
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)


@pytest.fixture(autouse=True)
def mock_plt_show(mocker):
    mocker.patch.object(plt, plt.show.__name__)


def run_example(example_name, expected_internal_energy, expected_error_message=None):
    if expected_error_message is not None:
        with pytest.raises(RuntimeError) as excinfo:
            importlib.import_module(example_name)
        assert str(excinfo.value) == expected_error_message
    else:
        importlib.import_module(example_name)
    assert_internal_energy(expected_internal_energy)


def assert_internal_energy(expected_internal_energy, tolerance=1e-5):
    db_handler = DataBaseHandler()
    for run_id, expected_uint in expected_internal_energy.items():
        uint = db_handler.scheme_tables.get_results(run_id)["uint"]
        assert abs(uint - expected_uint) < tolerance * abs(expected_uint)


def test_fixed_adr_qstls():
    expected_internal_energy = {
        1: -0.06915828119355959,
        2: -0.06915828119355959,
        3: -0.036272826904439975,
        4: -0.03545665197810883,
    }
    run_example("fixed_adr", expected_internal_energy)


def test_initial_guess_stls():
    expected_internal_energy = {1: -0.0696211823813233, 2: -0.06962120294563244}
    run_example("initial_guess_stls", expected_internal_energy)


def test_solve_quantum_schemes():
    expected_internal_energy = {1: -0.06915828119355959, 2: -0.07152825527327064}
    run_example("solve_quantum_schemes", expected_internal_energy)


def test_solve_qvsstls():
    expected_internal_energy = {8: -0.282671921632354}
    run_example("solve_qvsstls", expected_internal_energy)


def test_solve_rpa_and_esa():
    expected_internal_energy = {1: -0.09410276327343531, 2: -0.0700592697600796}
    run_example("solve_rpa_and_esa", expected_internal_energy)


def test_solve_stls():
    expected_internal_energy = {1: -0.0696212348180385}
    run_example("solve_stls", expected_internal_energy)


def test_solve_stls_iet():
    expected_internal_energy = {1: -0.07129865525802169, 2: -0.0716955500079742}
    run_example("solve_stls_iet", expected_internal_energy)


def test_solve_vsstls():
    expected_internal_energy = {19: -0.133311093051585}
    run_example("solve_vsstls", expected_internal_energy)
