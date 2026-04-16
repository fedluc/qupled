import subprocess
from pathlib import Path

from .uv import run_uv_tool


def docs():
    script_dir = Path(__file__).resolve().parent
    run_uv_tool(
        ["python", str(script_dir / "generate_diagrams.py")],
        extras=["docs"],
    )
    Path("docs", "_build", "doxygen").mkdir(parents=True, exist_ok=True)
    subprocess.run(["doxygen", "Doxyfile"], cwd="docs", check=True)
    run_uv_tool(
        ["sphinx-build", "-b", "html", "docs", str(Path("docs", "_build"))],
        extras=["docs"],
    )
