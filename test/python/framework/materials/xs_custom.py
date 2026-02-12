#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xs_custom.py: Regression test for user-provided custom XS loading.

Loads OpenMC HDPE cross sections, requests the "absorption" dataset by name
and compares the named xs values to the raw HDF5 dataset for a few indices.
"""

import os
import sys

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.xs import MultiGroupXS


def main():
    xs_file = "../../modules/linear_boltzmann_solvers/transport_steady/HDPE.h5"
    dataset = "set1"
    temperature = 294.0
    custom_xs = "absorption"
    tol = 1.0e-6
    expected = {
        0: 0.00434302,
        5: 0.00635482,
        9: 5.6538e-06,
        13: 5.05473e-06,
        17: 4.51771e-06,
        21: 4.1906e-06,
        25: 4.1563e-06,
        29: 5.00753e-06,
        33: 9.54322e-06,
        37: 1.80119e-05,
        41: 2.72368e-05,
        45: 4.56202e-05,
        49: 7.27354e-05,
        53: 0.000112619,
        57: 0.000150834,
        61: 0.000237293,
        65: 0.000412622,
        69: 0.000595737,
        73: 0.000700987,
        77: 0.000850868,
        81: 0.00113468,
        86: 0.00155062,
        91: 0.00216394,
        96: 0.00267416,
        101: 0.00303944,
        106: 0.00341594,
        110: 0.00367713,
        114: 0.00397723,
        118: 0.00415121,
        123: 0.00434302,
        128: 0.00454675,
        132: 0.00492324,
        136: 0.00604985,
        140: 0.00693216,
        144: 0.00786618,
        149: 0.0101541,
        154: 0.0133354,
        159: 0.0174687,
        164: 0.0263177,
        169: 0.0565075,
        171: 0.106801,
    }

    mgxs = MultiGroupXS()
    mgxs.LoadFromOpenMC(xs_file, dataset, temperature, [custom_xs])
    named_vals = mgxs.get_custom_xs(custom_xs)

    ok = True
    for idx, ref_val in expected.items():
        if idx >= len(named_vals):
            raise RuntimeError(f"Index {idx} out of range for {custom_xs}")
        diff = abs(float(named_vals[idx]) - ref_val)
        if diff > tol:
            ok = False
            break

    print(f"PASS {int(ok)}")


if __name__ == "__main__":
    main()
