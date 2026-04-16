import argparse

from .build import build
from .testing import test
from .install import install, install_dependencies
from .formatting import format_code
from .docs import docs
from .maintenance import clean, update_version


def run():
    parser = argparse.ArgumentParser(
        description="""A utility script for building, testing, formatting,
        and generating documentation for the qupled project."""
    )

    subparsers = parser.add_subparsers(dest="command", help="Sub-command to run")

    # Build command
    build_parser = subparsers.add_parser("build", help="Build the qupled package")
    build_parser.add_argument(
        "--use-mpi",
        action="store_true",
        help="Build with MPI support (default: False).",
    )
    build_parser.add_argument(
        "--native-only",
        action="store_true",
        help="Build only native code in C++ (default: False).",
    )

    # Test command
    test_parser = subparsers.add_parser("test", help="Run tests")
    test_parser.add_argument(
        "marker",
        nargs="?",
        choices=["unit", "native", "integration", "cpp"],
        default=None,
        help="Run only tests with this marker (default: run all tests).",
    )
    test_parser.add_argument(
        "--use-mpi",
        dest="use_mpi",
        action="store_true",
        default=False,
        help="Run cpp tests with MPI enabled.",
    )

    # Update version command
    version_parser = subparsers.add_parser(
        "update-version", help="Update package version"
    )
    version_parser.add_argument("build_version", help="The new version number.")

    # Other commands
    subparsers.add_parser("clean", help="Clean up build artifacts")
    subparsers.add_parser("docs", help="Generate documentation")
    subparsers.add_parser("format", help="Format the source code")
    subparsers.add_parser("install", help="Install the qupled package")
    subparsers.add_parser(
        "install-deps",
        help="Install system dependencies and sync Python dependencies with uv",
    )

    args = parser.parse_args()

    if args.command == "build":
        build(args.use_mpi, args.native_only)
    elif args.command == "clean":
        clean()
    elif args.command == "docs":
        docs()
    elif args.command == "format":
        format_code()
    elif args.command == "install":
        install()
    elif args.command == "test":
        test(args.marker, args.use_mpi)
    elif args.command == "install-deps":
        install_dependencies()
    elif args.command == "update-version":
        update_version(args.build_version)
    else:
        parser.print_help()


if __name__ == "__main__":
    run()
