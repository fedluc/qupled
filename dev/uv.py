import os
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
UV_EXTRA_FLAGS = ["--extra", "test", "--extra", "docs"]
UV_DEV_GROUP_FLAGS = ["--group", "dev"]
UV_RELEASE_GROUP_FLAGS = ["--group", "release"]


def run_uv(args, *, cwd=None, env=None):
    uv_path = shutil.which("uv")
    if uv_path is None:
        raise RuntimeError(
            "uv is required for Python environment management. "
            "Install it before running devtool."
        )
    subprocess.run([uv_path, *args], check=True, cwd=cwd, env=env)


def run_uv_tool(command, *, extras=None, groups=None, cwd=None, env=None):
    args = ["run"]
    args.extend(["--project", str(PROJECT_ROOT)])
    for extra in extras or []:
        args.extend(["--extra", extra])
    for group in groups or []:
        args.extend(["--group", group])
    args.extend(command)
    run_uv(args, cwd=cwd, env=env)


def sync_local_environment(*, skip_project_install=False):
    args = ["sync", *UV_DEV_GROUP_FLAGS, *UV_EXTRA_FLAGS]
    if skip_project_install:
        args.append("--no-install-project")
    run_uv(args, cwd=PROJECT_ROOT)


def sync_release_environment():
    run_uv(["sync", *UV_RELEASE_GROUP_FLAGS], cwd=PROJECT_ROOT)


def get_project_python():
    project_environment = os.environ.get("UV_PROJECT_ENVIRONMENT")
    if project_environment:
        environment_path = Path(project_environment)
    else:
        environment_path = PROJECT_ROOT / ".venv"

    if os.name == "nt":
        python_path = environment_path / "Scripts" / "python.exe"
    else:
        python_path = environment_path / "bin" / "python"
    if not python_path.exists():
        sync_local_environment()
    return str(python_path)


def print_uv_install_instructions():
    if sys.platform == "win32":
        print("Install uv from https://docs.astral.sh/uv/getting-started/installation/")
    else:
        print("Install uv with: curl -LsSf https://astral.sh/uv/install.sh | sh")
