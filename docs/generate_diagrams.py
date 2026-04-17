#!/usr/bin/env python3
"""Generate class inheritance diagrams (.dot files) for the documentation."""

import re
import sys
from pathlib import Path
from unittest.mock import MagicMock

ROOT = Path(__file__).parent.parent
INCLUDE_DIR = ROOT / "src" / "qupled" / "native" / "include"
OUTPUT_DIR = ROOT / "docs" / "_generated"

# ---------------------------------------------------------------------------
# Mock external C-extension and optional Python dependencies so that the
# scheme modules can be imported without a built package.
# ---------------------------------------------------------------------------
_MOCKS = [
    "qupled.native",
    "sqlalchemy",
    "sqlalchemy.orm",
    "sqlalchemy.dialects",
    "sqlalchemy.dialects.sqlite",
    "blosc2",
    "scipy",
    "scipy.integrate",
    "scipy.interpolate",
    "matplotlib",
    "matplotlib.pyplot",
]
for _mod in _MOCKS:
    if _mod not in sys.modules:
        sys.modules[_mod] = MagicMock()

sys.path.insert(0, str(ROOT / "src"))

from qupled.schemes import (  # noqa: E402
    esa,
    hf,
    qstls,
    qstlsiet,
    qvsstls,
    rpa,
    stls,
    stlsiet,
    vsstls,
)

# ---------------------------------------------------------------------------
# Python class lists
# ---------------------------------------------------------------------------
PYTHON_SOLVERS = [
    hf.Solver,
    rpa.Solver,
    esa.Solver,
    stls.Solver,
    stlsiet.Solver,
    vsstls.Solver,
    qstls.Solver,
    qstlsiet.Solver,
    qvsstls.Solver,
]
PYTHON_INPUTS = [
    hf.Input,
    rpa.Input,
    esa.Input,
    stls.Input,
    stlsiet.Input,
    vsstls.Input,
    qstls.Input,
    qstlsiet.Input,
    qvsstls.Input,
]

# ---------------------------------------------------------------------------
# C++ class lists (names only - parsed from headers)
# ---------------------------------------------------------------------------
CPP_SOLVERS = {
    "Logger",
    "HF",
    "Rpa",
    "ESA",
    "Stls",
    "StlsIet",
    "Qstls",
    "QstlsIet",
    "VSBase",
    "VSStls",
    "QVSStls",
}
CPP_INPUTS = {
    "Input",
    "IterationInput",
    "StlsInput",
    "StlsIetInput",
    "QstlsInput",
    "QstlsIetInput",
    "VSInput",
    "QuantumInput",
    "IetInput",
    "VSStlsInput",
    "QVSStlsInput",
}
CPP_MIXIN_CLASSES = {"VSInput", "QuantumInput", "IetInput"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_DOT_HEADER = """\
digraph {
    rankdir=BT
    node [shape=box, fontname="sans-serif", style=filled, fillcolor="#e8e8e8"]
    edge [arrowhead=empty]
"""
_DOT_FOOTER = "}\n"
_ROOT_COLOR = '"#c8daf5"'
_MIXIN_COLOR = '"#f5e8c8"'


def _short_name(cls):
    """Return 'module.ClassName' using the leaf module name."""
    mod = cls.__module__.split(".")[-1]
    return f"{mod}.{cls.__name__}"


def _build_python_dot(classes):
    class_set = set(classes)
    lines = [_DOT_HEADER]
    for cls in classes:
        name = _short_name(cls)
        if not any(b in class_set for b in cls.__bases__):
            lines.append(f'    "{name}" [fillcolor={_ROOT_COLOR}]')
    lines.append("")
    for cls in classes:
        name = _short_name(cls)
        for base in cls.__bases__:
            if base in class_set:
                lines.append(f'    "{name}" -> "{_short_name(base)}"')
    lines.append(_DOT_FOOTER)
    return "\n".join(lines)


def _parse_cpp_inheritance():
    """Return {class_name: [parent_names]} from all headers."""
    pattern = re.compile(
        r"class\s+(\w+)\s*:[^{]*?((?:(?:public|private|protected)\s+\w+\s*,?\s*)+)"
    )
    inheritance = {}
    for hpp in INCLUDE_DIR.rglob("*.hpp"):
        for match in pattern.finditer(hpp.read_text()):
            cls = match.group(1)
            parents = re.findall(
                r"(?:public|private|protected)\s+(\w+)", match.group(2)
            )
            inheritance[cls] = parents
    return inheritance


def _build_cpp_dot(target_classes, mixin_classes, inheritance):
    lines = [_DOT_HEADER]
    for cls in sorted(target_classes):
        parents_in_target = [p for p in inheritance.get(cls, []) if p in target_classes]
        if not parents_in_target:
            color = _MIXIN_COLOR if cls in mixin_classes else _ROOT_COLOR
            lines.append(f'    "{cls}" [fillcolor={color}]')
    for cls in mixin_classes:
        if cls in target_classes:
            lines.append(f'    "{cls}" [fillcolor={_MIXIN_COLOR}]')
    lines.append("")
    for cls in sorted(target_classes):
        for parent in inheritance.get(cls, []):
            if parent in target_classes:
                lines.append(f'    "{cls}" -> "{parent}"')
    lines.append(_DOT_FOOTER)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    (OUTPUT_DIR / "python_solvers.dot").write_text(_build_python_dot(PYTHON_SOLVERS))
    (OUTPUT_DIR / "python_inputs.dot").write_text(_build_python_dot(PYTHON_INPUTS))

    inheritance = _parse_cpp_inheritance()
    (OUTPUT_DIR / "cpp_solvers.dot").write_text(
        _build_cpp_dot(CPP_SOLVERS, set(), inheritance)
    )
    (OUTPUT_DIR / "cpp_inputs.dot").write_text(
        _build_cpp_dot(CPP_INPUTS, CPP_MIXIN_CLASSES, inheritance)
    )

    print("Diagrams generated in", OUTPUT_DIR)


if __name__ == "__main__":
    run()
