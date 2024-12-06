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
import sys

sys.path.insert(0, os.path.abspath(os.path.join("..", "src")))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.napoleon", "sphinx.ext.autodoc", "sphinx.ext.mathjax"]
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
