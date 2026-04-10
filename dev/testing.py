import os
import shutil
import subprocess
from pathlib import Path

from .build import NATIVE_TEST_BUILD_DIR, build_native_test_target
from .common import get_wheel_file


def run_tox(environment):
    tox_path = Path(".tox")
    if tox_path.exists():
        shutil.rmtree(tox_path)

    wheel_file = get_wheel_file()
    if wheel_file is not None:
        os.environ["WHEEL_FILE"] = wheel_file
        subprocess.run(["tox", "-e", environment], check=True)


def run_native_cpp_tests(use_mpi):
    build_native_test_target(use_mpi=use_mpi)
    subprocess.run(
        ["ctest", "--output-on-failure"], cwd=NATIVE_TEST_BUILD_DIR, check=True
    )


def test(marker, use_mpi):
    cpp_use_mpi = False if use_mpi is None else use_mpi

    if marker is None:
        run_tox("test")
        run_native_cpp_tests(use_mpi=cpp_use_mpi)
    elif marker == "cpp":
        run_native_cpp_tests(use_mpi=cpp_use_mpi)
    else:
        run_tox(marker)
