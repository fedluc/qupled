import os
import shutil
import subprocess
from pathlib import Path

from .common import get_wheel_file
from .uv import (
    get_project_python,
    print_uv_install_instructions,
    sync_local_environment,
)


def install():
    wheel_file = get_wheel_file()
    if wheel_file is not None:
        sync_local_environment(skip_project_install=True)
        python_path = get_project_python()
        subprocess.run(
            ["uv", "pip", "uninstall", "--python", python_path, "qupled"],
            check=False,
        )
        subprocess.run(
            ["uv", "pip", "install", "--python", python_path, wheel_file], check=True
        )


def install_dependencies():
    print("Installing dependencies...")
    script_dir = Path(__file__).resolve().parent
    if os.name == "posix":
        if shutil.which("apt-get"):
            _install_with_apt(script_dir / "requirements-apt.txt")
        elif shutil.which("brew"):
            _install_with_brew(script_dir / "requirements-brew.txt")
        else:
            print("Unsupported package manager. Please install dependencies manually.")
    else:
        print("Unsupported operating system. Please install dependencies manually.")
    try:
        sync_local_environment(skip_project_install=True)
    except RuntimeError as exc:
        print(exc)
        print_uv_install_instructions()
        raise SystemExit(1) from exc


def _install_with_apt(apt_requirements):
    subprocess.run(["sudo", "apt-get", "update"], check=True)
    with apt_requirements.open("r") as apt_file:
        subprocess.run(
            ["xargs", "sudo", "apt-get", "install", "-y"],
            stdin=apt_file,
            check=True,
        )


def _install_with_brew(brew_requirements):
    subprocess.run(["brew", "update"], check=True)
    subprocess.run(["brew", "bundle", f"--file={brew_requirements}"], check=True)
