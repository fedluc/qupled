import subprocess
import platform

def install_dependencies():
    system = platform.system()
    try:
        if system == "Linux":
            print("Updating package list...")
            subprocess.check_call(["apt-get", "update"])
            print("Installing dependencies...")
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
            print("Installing dependencies using Homebrew...")
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
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        print(f"Command: {e.cmd}")
        if e.output:
            print(f"Output: {e.output.decode('utf-8')}")
        raise
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise

if __name__ == "__main__":
    install_dependencies()
