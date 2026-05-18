#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import runpy
from pathlib import Path

script_globals = dict(globals())
script_globals["groupset_mode"] = "one_group_per_groupset"

runpy.run_path(
    str(Path(__file__).with_name("c5g7.py")),
    run_name="__main__",
    init_globals=script_globals,
)
