import importlib
import os
import sys

import json
import matplotlib.pyplot as plt
import pytest

from qupled.database import DataBaseHandler


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
    database_name = DataBaseHandler.DEFAULT_DATABASE_NAME
    if os.path.exists(database_name):
        os.remove(database_name)


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


def assert_internal_energy(expected_internal_energy, tolerance=1e-10):
    db_handler = DataBaseHandler()
    for run_id, expected_uint in expected_internal_energy.items():
        uint = db_handler.get_results(run_id)["uint"]
        assert abs(uint - expected_uint) < tolerance


def test_fixed_adr_qstls():
    expected_internal_energy = {
        1: -0.06915834235189136,
        2: -0.06915834235189136,
        3: -0.03627284914113849,
        4: -0.035456654705259216,
    }
    run_example("fixed_adr", expected_internal_energy)


def test_initial_guess_stls():
    expected_internal_energy = {1: -0.06962128794619024, 2: -0.06962128710276651}
    run_example("initial_guess_stls", expected_internal_energy)


def test_solve_quantum_schemes():
    expected_internal_energy = {1: -0.06915834235189136, 2: -0.07152838527583957}
    run_example("solve_quantum_schemes", expected_internal_energy)


def test_solve_qvsstls():
    expected_internal_energy = {
        1: -3.021373777177104,
        2: -1.0023714941672779,
        3: -0.6439252234345055,
        4: -0.4840844703408011,
        5: -0.39133370228089476,
        6: -0.3299820675074039,
        7: -0.28606359879189935,
        8: -0.274026991680492,
    }
    run_example("solve_qvsstls", expected_internal_energy)


def test_solve_rpa_and_esa():
    expected_internal_energy = {1: -0.09410276327343531, 2: -0.07005926976007962}
    run_example("solve_rpa_and_esa", expected_internal_energy)


def test_solve_stls():
    expected_internal_energy = {1: -0.06962128635463076}
    run_example("solve_stls", expected_internal_energy)


def test_solve_stls_iet():
    expected_internal_energy = {1: -0.07129872557238329, 2: -0.07169561914217}
    run_example("solve_stls_iet", expected_internal_energy)


def test_solve_vsstls():
    expected_internal_energy = {
        1: -2.9940880309281104,
        2: -0.9993088491331527,
        3: -0.6469894140641651,
        4: -0.490134772126587,
        5: -0.3990831438584432,
        6: -0.33876127750338797,
        7: -0.29547372271661176,
        8: -0.28358922791731095,
        9: -0.262698528855232,
        10: -0.23689815309397158,
        11: -0.21599479260943763,
        12: -0.19866597531856445,
        13: -0.18402529485178737,
        14: -0.17149705471535312,
        15: -0.16062058250083427,
        16: -0.1511518189162366,
        17: -0.14288119657620452,
        18: -0.13553308850917165,
        19: -0.13334619793591626,
    }
    run_example("solve_vsstls", expected_internal_energy)
