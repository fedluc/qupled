import importlib
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from qupled.database.database_handler import DATABASE_DIRECTORY


@pytest.fixture(autouse=True)
def run_before_each_test():
    examples_dir = str(Path(__file__).resolve().parents[2] / "examples" / "docs")
    if examples_dir not in sys.path:
        sys.path.insert(0, examples_dir)
    yield
    if examples_dir in sys.path:
        sys.path.remove(examples_dir)


@pytest.fixture(autouse=True)
def run_after_each_test():
    yield
    output_dir = DATABASE_DIRECTORY
    if Path(output_dir).exists():
        shutil.rmtree(output_dir)


@pytest.fixture(autouse=True)
def mock_plt_show(mocker):
    mocker.patch.object(plt, plt.show.__name__)


def run_example(example_name):
    importlib.import_module(example_name)


@pytest.mark.integration
def test_finite_size_correction():
    run_example("finite_size_correction")


@pytest.mark.integration
def test_fixed_adr_qstls():
    run_example("fixed_adr")


@pytest.mark.integration
def test_initial_guess_stls():
    run_example("initial_guess_stls")


@pytest.mark.integration
def test_post_processing():
    run_example("post_processing")


@pytest.mark.integration
def test_solve_quantum_schemes():
    run_example("solve_quantum_schemes")


@pytest.mark.integration
def test_solve_rpa_and_esa():
    run_example("solve_rpa_and_esa")


@pytest.mark.integration
def test_solve_stls():
    run_example("solve_stls")


@pytest.mark.integration
def test_solve_vsstls():
    run_example("solve_vsstls")
