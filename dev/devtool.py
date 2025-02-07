import os
import argparse
import subprocess
import shutil
from pathlib import Path


def build(nompi, native_only):
    # Build without MPI
    if nompi:
        os.environ["USE_MPI"] = "OFF"
    # Set environment variable for OpenMP on macOS
    if os.name == "posix" and shutil.which("brew"):
        brew_prefix = subprocess.run(
            ["brew", "--prefix"], capture_output=True, text=True
        ).stdout.strip()
        os.environ["OpenMP_ROOT"] = str(Path(brew_prefix, "opt", "libomp"))
    if native_only:
        build_folder = "dist-native-only"
        if not os.path.exists(build_folder):
            os.makedirs(build_folder)
        os.chdir(build_folder)
        subprocess.run(["cmake", "../src/qupled/native/src"], check=True)
        subprocess.run(["cmake", "--build", "."], check=True)
    else:
        subprocess.run(["python3", "-m", "build"], check=True)
    print("Build completed.")


def get_wheel_file():
    wheel_file = list(Path().rglob("qupled*.whl"))
    if not wheel_file:
        print("No .whl files found. Ensure the package is built first.")
        return None
    else:
        return str(wheel_file[0])


def run_tox(environment):
    tox_path = Path(".tox")
    if tox_path.exists():
        shutil.rmtree(tox_path)
    wheel_file = get_wheel_file()
    if wheel_file is not None:
        os.environ["WHEEL_FILE"] = wheel_file
        subprocess.run(["tox", "-e", environment], check=True)


def test():
    run_tox("test")


def examples():
    run_tox("examples")


def format_code():
    subprocess.run(["black", "."], check=True)
    native_files_folder = Path("src", "qupled", "native")
    cpp_files = list(native_files_folder.rglob("*.cpp"))
    hpp_files = list(native_files_folder.rglob("*.hpp"))
    for f in cpp_files + hpp_files:
        subprocess.run(["clang-format", "--style=file", "-i", str(f)], check=True)


def docs():
    subprocess.run(["sphinx-build", "-b", "html", "docs", str(Path("docs", "_build"))])


def clean():
    folders_to_clean = [
        Path("dist"),
        Path("dist-native-only"),
        Path("src", "qupled.egg-info"),
        Path("docs", "_build"),
    ]
    for folder in folders_to_clean:
        if folder.exists():
            print(f"Removing folder: {folder}")
            shutil.rmtree(folder)


def install():
    wheel_file = get_wheel_file()
    if wheel_file is not None:
        subprocess.run(["pip", "install", "--force-reinstall", wheel_file], check=True)


def install_dependencies():
    print("Installing dependencies...")
    if os.name == "posix":
        if shutil.which("apt-get"):
            subprocess.run(["sudo", "apt-get", "update"], check=True)
            subprocess.run(
                [
                    "sudo",
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
                ],
                check=True,
            )
        elif shutil.which("brew"):
            subprocess.run(["brew", "update"], check=True)
            subprocess.run(
                [
                    "brew",
                    "install",
                    "cmake",
                    "gsl",
                    "libomp",
                    "openmpi",
                    "fmt",
                    "boost-python3",
                ],
                check=True,
            )
        else:
            print("Unsupported package manager. Please install dependencies manually.")
    else:
        print("Unsupported operating system. Please install dependencies manually.")


def update_version(build_version):
    pyproject_file = Path("pyproject.toml")
    if not pyproject_file.exists():
        return
    with pyproject_file.open("r") as file:
        content = file.readlines()
    with pyproject_file.open("w") as file:
        for line in content:
            if line.startswith("version = "):
                file.write(f'version = "{build_version}"')
                file.write("\n")
            else:
                file.write(line)


def run():
    parser = argparse.ArgumentParser(
        description="""A utility script for building, testing, formatting,
        and generating documentation for the qupled project."""
    )

    subparsers = parser.add_subparsers(dest="command", help="Sub-command to run")

    # Build command
    build_parser = subparsers.add_parser("build", help="Build the qupled package")
    build_parser.add_argument(
        "--nompi",
        action="store_true",
        help="Build without MPI support (default: False).",
    )
    build_parser.add_argument(
        "--native-only",
        action="store_true",
        help="Build only native code in C++ (default: False).",
    )

    # Update version command
    version_parser = subparsers.add_parser(
        "update-version", help="Update package version"
    )
    version_parser.add_argument("build_version", help="The new version number.")

    # Other commands
    subparsers.add_parser("clean", help="Clean up build artifacts")
    subparsers.add_parser("docs", help="Generate documentation")
    subparsers.add_parser("examples", help="Run tests for the examples")
    subparsers.add_parser("format", help="Format the source code")
    subparsers.add_parser("install", help="Install the qupled package")
    subparsers.add_parser("install-deps", help="Install system dependencies")
    subparsers.add_parser("test", help="Run tests")

    args = parser.parse_args()

    if args.command == "build":
        build(args.nompi, args.native_only)
    elif args.command == "clean":
        clean()
    elif args.command == "docs":
        docs()
    elif args.command == "examples":
        examples()
    elif args.command == "format":
        format_code()
    elif args.command == "install":
        install()
    elif args.command == "test":
        test()
    elif args.command == "install-deps":
        install_dependencies()
    elif args.command == "update-version":
        update_version(args.build_version)
    else:
        parser.print_help()


if __name__ == "__main__":
    run()
