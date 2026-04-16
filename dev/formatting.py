import subprocess
from pathlib import Path

from .uv import run_uv_tool


def format_code():
    run_uv_tool(["black", "."], groups=["dev"])
    native_files_folder = Path("src", "qupled", "native")
    cpp_files = list(native_files_folder.rglob("*.cpp"))
    hpp_files = list(native_files_folder.rglob("*.hpp"))
    for file_path in cpp_files + hpp_files:
        subprocess.run(
            ["clang-format", "--style=file", "-i", str(file_path)], check=True
        )
