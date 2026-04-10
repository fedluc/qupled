import subprocess
from pathlib import Path


def docs():
    script_dir = Path(__file__).resolve().parent
    subprocess.run(["python3", str(script_dir / "generate_diagrams.py")], check=True)
    Path("docs", "_build", "doxygen").mkdir(parents=True, exist_ok=True)
    subprocess.run(["doxygen", "Doxyfile"], cwd="docs", check=True)
    subprocess.run(["sphinx-build", "-b", "html", "docs", str(Path("docs", "_build"))])
