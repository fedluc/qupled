import os
import argparse
import subprocess
import shutil
from pathlib import Path


def build(nompi):
    if nompi:
        os.environ["USE_MPI"] = "OFF"
    subprocess.run(["python3", "-m", "build"], check=True)
    print("Build completed.")


def get_wheel_file():
    wheel_file = list(Path(".").rglob("qupled*.whl"))
    if not wheel_file:
        print("No .whl files found. Ensure the package is built first.")
        return None
    else:
        return str(wheel_file[0])


def run_tox(environment):
    if os.path.exists(".tox"):
        shutil.rmtree(".tox")
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
    cpp_files = list(Path("qupled").rglob("*.cpp"))
    hpp_files = list(Path("qupled").rglob("*.hpp"))
    for f in cpp_files + hpp_files:
        subprocess.run(["clang-format", "--style=file", "-i", str(f)], check=True)


def docs():
    subprocess.run(
        ["sphinx-build", "-b", "html", "docs", os.path.join("docs", "_build")]
    )


def clean():
    folders_to_clean = [
        "dist",
        os.path.join("src", "qupled.egg-info"),
        os.path.join("docs", "_build"),
    ]
    for folder in folders_to_clean:
        if os.path.exists(folder):
            print(f"Removing folder: {folder}")
            shutil.rmtree(folder)


def install():
    wheel_file = get_wheel_file()
    if wheel_file is not None:
        subprocess.run(["pip", "install", "--force-reinstall", wheel_file], check=True)


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

    # Other commands
    subparsers.add_parser("clean", help="Clean up build artifacts")
    subparsers.add_parser("docs", help="Generate documentation")
    subparsers.add_parser("examples", help="Run tests for the examples")
    subparsers.add_parser("format", help="Format the source code")
    subparsers.add_parser("install", help="Install the qupled package")
    subparsers.add_parser("test", help="Run tests")

    args = parser.parse_args()

    if args.command == "build":
        build(args.nompi)
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
    else:
        parser.print_help()
