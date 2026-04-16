import subprocess

from .build import NATIVE_TEST_BUILD_DIR, build_native_test_target
from .uv import run_uv_tool


def run_pytest(arguments):
    run_uv_tool(
        ["pytest", *arguments],
        extras=["test"],
    )


def run_native_cpp_tests(use_mpi):
    build_native_test_target(use_mpi=use_mpi)
    subprocess.run(
        ["ctest", "--output-on-failure"], cwd=NATIVE_TEST_BUILD_DIR, check=True
    )


def test(marker, use_mpi):
    cpp_use_mpi = False if use_mpi is None else use_mpi

    if marker is None:
        run_pytest(["tests", "-v"])
        run_native_cpp_tests(use_mpi=cpp_use_mpi)
    elif marker == "unit":
        run_pytest(["-m", "unit", "tests", "-v"])
    elif marker == "native":
        run_pytest(["-m", "native", "tests", "-v"])
    elif marker == "integration":
        run_pytest(["-m", "integration", "tests", "-v"])
    elif marker == "cpp":
        run_native_cpp_tests(use_mpi=cpp_use_mpi)
