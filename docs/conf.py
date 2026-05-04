# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "qupled"
copyright = "2023, Federico Lucco Castello"
author = "Federico Lucco Castello"

# -- Path to modules ---------------------------------------------------------
import os
import subprocess
import sys

sys.path.insert(0, os.path.abspath(os.path.join("..", "src")))

# -- Run Doxygen before Breathe reads the XML --------------------------------
_docs_dir = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(_docs_dir, "_build", "doxygen"), exist_ok=True)
subprocess.run(["doxygen", "Doxyfile"], cwd=_docs_dir, check=True)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.graphviz",
    "breathe",
    "sphinxcontrib.mermaid",
]

breathe_projects = {"qupled": "_build/doxygen/xml"}
breathe_default_project = "qupled"
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
autodoc_member_order = "bysource"
autodoc_mock_imports = ["qupled.native"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = "sphinx_rtd_theme"
html_theme_options = {"collapse_navigation": False, "navigation_depth": 3}
html_static_path = ["_static"]
html_css_files = ["css/rdt_theme_python_properties.css"]
html_show_sourcelink = False
