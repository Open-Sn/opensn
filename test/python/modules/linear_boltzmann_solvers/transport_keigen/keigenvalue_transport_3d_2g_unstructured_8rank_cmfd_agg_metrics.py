#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8-rank unstructured CMFD regression with verbose coarse-mesh diagnostics enabled.

This regression reuses the 8-rank extruded-unstructured CMFD k-eigenvalue
benchmark and enables verbose CMFD diagnostics inline so the regression harness
can inspect the rank-local coarse-mesh aggregation. It verifies the benchmark
k-eigenvalue, a bounded accelerated sweep count, and the reported aggregation
ratio for the default unstructured mesh configuration.
"""

import os
import runpy

cmfd_verbose = True

runpy.run_path(
    os.path.join(
        os.path.dirname(__file__),
        "keigenvalue_transport_3d_2g_unstructured_8rank_cmfd_agg.py",
    ),
    init_globals=globals(),
    run_name="__main__",
)
