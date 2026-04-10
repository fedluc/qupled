import os
import shutil
import subprocess
from pathlib import Path

NATIVE_BUILD_DIR = Path("dist-native-only")
NATIVE_SOURCE_DIR = Path("src", "qupled", "native", "src")


def configure_openmp_root_for_macos():
    if os.name != "posix" or not shutil.which("brew"):
        return
    if "OpenMP_ROOT" in os.environ:
        return

    brew_prefix = subprocess.run(
        ["brew", "--prefix"], capture_output=True, text=True, check=True
    ).stdout.strip()
    os.environ["OpenMP_ROOT"] = str(Path(brew_prefix, "opt", "libomp"))


def build(use_mpi, native_only, native_tests):
    # Build with MPI
    if use_mpi:
        os.environ["USE_MPI"] = "ON"
    # Build native C++ unit tests
    if native_tests:
        os.environ["BUILD_NATIVE_TESTS"] = "ON"
    # Set environment variable for OpenMP on macOS
    configure_openmp_root_for_macos()
    if native_only:
        build_native(native_tests=native_tests)
    else:
        subprocess.run(["python3", "-m", "build"], check=True)

    print("Build completed.")


def build_native(native_tests):
    configure_openmp_root_for_macos()
    NATIVE_BUILD_DIR.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "cmake",
            str(NATIVE_SOURCE_DIR.resolve()),
            f"-DBUILD_NATIVE_TESTS={'ON' if native_tests else 'OFF'}",
        ],
        cwd=NATIVE_BUILD_DIR,
        check=True,
    )
    subprocess.run(["cmake", "--build", "."], cwd=NATIVE_BUILD_DIR, check=True)
