#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xs_custom_combine_scale.py: Regression test for custom XS combine and scale behavior.

Loads a named 1D XS from OpenMC data, verifies that scaling affects the custom XS values,
and verifies that combining two cross sections preserves and density-weights the custom XS.
"""

import os
import sys

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.xs import MultiGroupXS

xs_file = "../../modules/linear_boltzmann_solvers/transport_steady/HDPE.h5"
dataset = "set1"
temperature = 294.0
custom_xs = "absorption"
probe_index = 81
tol = 1.0e-6
expected = {
    "scaled": 0.0028367,
    "combined": 0.00397138,
}

xs_1 = MultiGroupXS()
xs_1.LoadFromOpenMC(xs_file, dataset, temperature, [custom_xs])
xs_1.Scale(2.5)
scaled_value = float(xs_1.get_custom_xs(custom_xs)[probe_index])

xs_2 = MultiGroupXS()
xs_2.LoadFromOpenMC(xs_file, dataset, temperature, [custom_xs])
xs_3 = MultiGroupXS()
xs_3.LoadFromOpenMC(xs_file, dataset, temperature, [custom_xs])

d1 = 0.5
d2 = 3.0
xs_combined = MultiGroupXS.Combine([(xs_2, d1), (xs_3, d2)])
combined_value = float(xs_combined.get_custom_xs(custom_xs)[probe_index])

ok = True
if abs(scaled_value - expected["scaled"]) > tol:
    ok = False
if abs(combined_value - expected["combined"]) > tol:
    ok = False

print(f"PASS {int(ok)}")
