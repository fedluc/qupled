import importlib
import os
import sys

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
    expected_internal_energy = {1: -0.06962131837263942, 2: -0.06962131103269621}
    run_example("initial_guess_stls", expected_internal_energy)


def test_solve_quantum_schemes():
    expected_internal_energy = {1: -0.06915834235189136, 2: -0.07152833875922872}
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
    expected_internal_energy = {1: -0.06962129065816987}
    run_example("solve_stls", expected_internal_energy)


def test_solve_stls_iet():
    expected_internal_energy = {1: -0.07129873381329867, 2: -0.0716956405760532}
    run_example("solve_stls_iet", expected_internal_energy)


def test_solve_vsstls():
    expected_internal_energy = {
        1: -2.994095202019811,
        2: -0.9993104835032114,
        3: -0.6469903999058134,
        4: -0.490135344559391,
        5: -0.3990834406333417,
        6: -0.3387613792997866,
        7: -0.2954736780653446,
        8: -0.2835891671524341,
        9: -0.26269815165830607,
        10: -0.23689855826562306,
        11: -0.21599398948537632,
        12: -0.1986644854596961,
        13: -0.18402552349746606,
        14: -0.17149382511359426,
        15: -0.1606191478690345,
        16: -0.15114500705443706,
        17: -0.14287185591472557,
        18: -0.1355613009089871,
        19: -0.13329527133219166,
    }
    run_example("solve_vsstls", expected_internal_energy)
