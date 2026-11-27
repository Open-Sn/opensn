# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import logging
import re
import sys
from pathlib import Path

# -- Import package ----------------------------------------------------------

project_dir = Path(__file__).resolve().parent.parent.parent
binary_dir = project_dir / "build"

# Prefer the newly built package (build/pyopensn) over the project directory
# or any site-packages install. Remove-then-insert so the final order is deterministic
# even when one of these paths is already on sys.path, e.g. via PYTHONPATH.
for path in (str(project_dir), str(binary_dir)):
    if path in sys.path:
        sys.path.remove(path)
    sys.path.insert(0, path)

import pyopensn  # noqa: E402

# The coupling module exists only in builds configured with OPENSN_WITH_PRECICE=ON. Without it,
# document a placeholder so the API index still builds.
if not hasattr(pyopensn, "coupling"):
    import types

    _coupling = types.ModuleType("pyopensn.coupling", "Multiphysics coupling module.")

    class PreciceAdapter:
        """One-way OpenSn preCICE volume-coupling participant.

        Not available in this build. Configure OpenSn with ``OPENSN_WITH_PRECICE=ON`` to build
        it; see the preCICE coupling section of the user guide.
        """

    PreciceAdapter.__module__ = "pyopensn.coupling"
    _coupling.PreciceAdapter = PreciceAdapter
    pyopensn.coupling = _coupling
    sys.modules["pyopensn.coupling"] = _coupling

# -- Project information -----------------------------------------------------

project = "OpenSn"
copyright = "2023-present, OpenSn team"
author = "OpenSn team"
release = (project_dir / "VERSION.txt").read_text().strip()
if pyopensn.__version__ != release:
    raise RuntimeError(
        f"The pyopensn build is version {pyopensn.__version__}, but VERSION.txt contains {release}"
    )

# Dependency floors come from the file CMake reads, as |min_petsc| and so on.
_floors_text = (project_dir / "cmake" / "OpenSnDependencyVersions.cmake").read_text()
_floor_lines = re.findall(r"^\s*set\(OPENSN_MIN_", _floors_text, re.MULTILINE)
_floors = re.findall(r"^set\(OPENSN_MIN_(\w+)_VERSION ([0-9.]+)\)", _floors_text, re.MULTILINE)
if not _floors or len(_floors) != len(_floor_lines):
    raise RuntimeError("Could not read every floor in cmake/OpenSnDependencyVersions.cmake")
rst_prolog = "\n".join(
    f".. |min_{name.lower()}| replace:: {version}" for name, version in _floors
)

# -- General configuration ---------------------------------------------------

master_doc = "index"
templates_path = ["_templates"]
source_suffix = [".rst"]
exclude_patterns = ["tutorials/_templates/*"]

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
