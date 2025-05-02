import glob
import importlib
import os
import sys

import matplotlib.pyplot as plt
import pytest

from qupled.database import DataBaseHandler


@pytest.fixture(scope="session", autouse=True)
def setup_examples_dir():
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


def run_example(example_name, expected_error_message=None):
    if expected_error_message is not None:
        with pytest.raises(RuntimeError) as excinfo:
            importlib.import_module(example_name)
        assert str(excinfo.value) == expected_error_message
    else:
        importlib.import_module(example_name)


def test_fixed_adr_qstls():
    run_example("fixed_adr")


def test_initial_guess_stls():
    run_example("initial_guess_stls")


def test_solve_quantum_schemes():
    run_example("solve_quantum_schemes")


def test_solve_qvs_stls():
    run_example("solve_qvsstls")


def test_solve_rpa_and_esa():
    run_example("solve_rpa_and_esa")


def test_solve_stls():
    run_example("solve_stls")


def test_solve_stls_iet():
    run_example("solve_stls_iet")


def test_solve_vs_stls():
    run_example("solve_vsstls")
