# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import logging
import sys
from pathlib import Path

# -- Import package ----------------------------------------------------------

project_dir = Path(__file__).resolve().parent.parent.parent
binary_dir = project_dir / "build"

# Prefer the in-tree package over any stale site-packages install.
for path in (str(project_dir), str(binary_dir)):
    if path not in sys.path:
        sys.path.insert(0, path)

import pyopensn  # noqa: E402

# -- Project information -----------------------------------------------------

project = "OpenSn"
copyright = "2023-present, OpenSn team"
author = "OpenSn team"
release = pyopensn.__version__

# -- General configuration ---------------------------------------------------

master_doc = "index"
templates_path = ["_templates"]
source_suffix = [".rst"]

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinxcontrib.bibtex",
    "sphinx_design",
    "nbsphinx",
    "breathe",
    "exhale",
]

nbsphinx_execute = "never"

breathe_projects = {
    "opensn": str(Path(__file__).resolve().parent / "capi" / "xml")
}
breathe_default_project = "opensn"

exhale_args = {
    "containmentFolder": str(Path(__file__).resolve().parent / "capi" / "rst"),
    "rootFileName": "index.rst",
    "rootFileTitle": "C++ API",
    "doxygenStripFromPath": str(Path(__file__).resolve().parent.parent.parent),
    "exhaleExecutesDoxygen": False,
    "createTreeView": True,
    "fullToctreeMaxDepth": 1,
    "unabridgedOrphanKinds": {"dir", "file", "page", "class", "enum", "namespace"},
    "fullApiSubSectionTitle": "Extras"
}

cpp_id_attributes = ['__inline_host_dev__']

# -- Options for BibTeX ------------------------------------------------------

bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "label"

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_show_sourcelink = False
html_css_files = ["pyopensn.css"]

# -- Options for LaTeX output -------------------------------------------------

latex_elements = {
    "maxlistdepth": "10"
}


class _CppWarningFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        message = record.getMessage()
        suppressed_fragments = (
            "Duplicate C++ declaration",
            "Invalid C++ declaration",
            "Error when parsing function declaration.",
            "Candidate function could not be parsed.",
        )
        return not any(fragment in message for fragment in suppressed_fragments)


def setup(app):
    logging.getLogger("sphinx").addFilter(_CppWarningFilter())
