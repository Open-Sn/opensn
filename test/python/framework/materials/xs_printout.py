#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read a cross section file, create a scaled copy, and print out some components to the console.
"""

import pprint
import logging


def dump(o):
    return pprint.pformat(o)

# Create cross sections
# from pyopensn.xs import MultiGroupXS


my_xs = {}

# Load cross section data for fuel
fuel_xs = MultiGroupXS()
fuel_xs.LoadFromOpenMC(
    "../../modules/linear_boltzmann_solvers/transport_keigen/u235_172g.h5",
    "u235",
    294.0
)
my_xs["fuel"] = fuel_xs

# Before scaling
chi_before = my_xs["fuel"].chi[0]
sigt_before = my_xs["fuel"].sigma_t[0]

# Create a scaled copy
my_xs["fuel"].Scale(2.0)

# After scaling
chi_after = my_xs["fuel"].chi[0]
sigt_after = my_xs["fuel"].sigma_t[0]

# Output
print(f"chi[0] before: {chi_before}")
print(f"chi[0] after : {chi_after}")
print(f"sigt[0] before: {sigt_before}")
print(f"sigt[0] after : {sigt_after}")
