import os
import shutil
import subprocess
from pathlib import Path

from .uv import run_uv_tool

NATIVE_BUILD_DIR = Path("dist-native")
NATIVE_TEST_BUILD_DIR = Path("dist-native-tests")
NATIVE_SOURCE_DIR = Path("src", "qupled", "native", "src")


def build(use_mpi, native_only):
    _configure_openmp_root_for_macos()
    if native_only:
        build_native(native_tests=False, use_mpi=use_mpi)
    else:
        _build_python_package(use_mpi=use_mpi)
    print("Build completed.")


def _to_cmake_bool(value):
    return "ON" if value else "OFF"


def _configure_openmp_root_for_macos():
    if os.name != "posix" or not shutil.which("brew"):
        return
    if "OpenMP_ROOT" in os.environ:
        return
    brew_prefix = subprocess.run(
        ["brew", "--prefix"], capture_output=True, text=True, check=True
    ).stdout.strip()
    os.environ["OpenMP_ROOT"] = str(Path(brew_prefix, "opt", "libomp"))


def _native_cmake_configure(native_tests, use_mpi, build_dir):
    subprocess.run(
        [
            "cmake",
            str(NATIVE_SOURCE_DIR.resolve()),
            f"-DBUILD_NATIVE_TESTS={_to_cmake_bool(native_tests)}",
            f"-DUSE_MPI={_to_cmake_bool(use_mpi)}",
        ],
        cwd=build_dir,
        check=True,
    )


def _native_cmake_build(build_dir, target=None):
    command = ["cmake", "--build", ".", "--parallel"]
    if target is not None:
        command.extend(["--target", target])
    subprocess.run(command, cwd=build_dir, check=True)


def _build_python_package(use_mpi):
    build_env = os.environ.copy()
    build_env["USE_MPI"] = _to_cmake_bool(use_mpi)
    build_env["BUILD_NATIVE_TESTS"] = "OFF"

    run_uv_tool(
        ["python", "-m", "build"],
        groups=["dev"],
        env=build_env,
    )


def build_native(native_tests, use_mpi):
    _configure_openmp_root_for_macos()
    NATIVE_BUILD_DIR.mkdir(parents=True, exist_ok=True)
    _native_cmake_configure(
        native_tests=native_tests, use_mpi=use_mpi, build_dir=NATIVE_BUILD_DIR
    )
    _native_cmake_build(build_dir=NATIVE_BUILD_DIR)


def build_native_test_target(use_mpi=True):
    _configure_openmp_root_for_macos()
    NATIVE_TEST_BUILD_DIR.mkdir(parents=True, exist_ok=True)
    _native_cmake_configure(
        native_tests=True, use_mpi=use_mpi, build_dir=NATIVE_TEST_BUILD_DIR
    )
    _native_cmake_build(build_dir=NATIVE_TEST_BUILD_DIR, target="native_tests")
