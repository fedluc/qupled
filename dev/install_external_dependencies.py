import os
import platform
import subprocess


def install_dependencies():
    system = platform.system()
    if system == "Linux":
        # subprocess.check_call(["apt-get", "update"])
        subprocess.check_call(
            [
                "apt-get",
                "install",
                "-y",
                "cmake",
                "libboost-all-dev",
                "libopenmpi-dev",
                "libgsl-dev",
                "libomp-dev",
                "libfmt-dev",
                "python3-dev",
            ]
        )
    elif system == "Darwin":
        subprocess.check_call(
            [
                "brew",
                "install",
                "cmake",
                "gsl",
                "libomp",
                "openmpi",
                "fmt",
                "boost-python3",
            ]
        )
    else:
        raise RuntimeError(f"Unsupported platform: {system}")


if __name__ == "__main__":
    install_dependencies()
