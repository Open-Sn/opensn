#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 2G heterogeneous qblock CMFD regression with all fine energy groups collapsed
into one coarse CMFD group.

This regression reuses the standard spatially aggregated 3D 2-group qblock
CMFD k-eigenvalue problem, but sets group_aggregation_size=2 so both transport
groups collapse into one coarse CMFD energy group. It checks that energy-group
collapse is wired into coarse-mesh setup, that the verbose coarse-mesh
diagnostic reports one coarse energy group, and that the accelerated solve
preserves the expected k-eigenvalue.
"""

import os
import runpy

cmfd_group_aggregation_size = 2
cmfd_verbose = True

runpy.run_path(
    os.path.join(os.path.dirname(__file__), "keigenvalue_transport_3d_2g_qblock_cmfd_agg.py"),
    init_globals=globals(),
    run_name="__main__",
)
