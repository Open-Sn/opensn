# Tutorial authoring template

Copy `tutorial_template.ipynb` into the appropriate subject directory and
replace every bracketed placeholder. Keep a tutorial focused on one learning
goal and prefer data already maintained by the regression suite.

Each executable notebook should:

- state its audience, prerequisites, learning objectives, and expected output;
- use paths relative to the notebook rather than user-specific absolute paths;
- finish with a stable, quantitative verification and a machine-readable output;
- have a neighboring `tests.json` entry; and
- be included in the nearest `index.rst` toctree.