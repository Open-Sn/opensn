#!/usr/bin/env python3
"""Regression test for 3D interpolation of the OpenMC 2a fuel XS data."""

import numpy as np

if "opensn_console" not in globals():
    from pyopensn.xs import CartesianGrid, LinearInterpolator, MultiGroupXS, XSType


# The expected values were independently calculated from the eight OpenMC
# corners with trilinear weights at STATE_POINT.
GRID_DATA = [[600.0, 1300.0], [900.0, 1200.0], [500.0, 700.0]]
STATE_POINT = np.array([990.0, 1040.0, 620.0], dtype=float)
CORNER_FILES = {
    (0, 0, 0): "../../../assets/xs/mgxs_2a_Bo6Tf9Tm5.h5",
    (0, 0, 1): "../../../assets/xs/mgxs_2a_Bo6Tf9Tm7.h5",
    (0, 1, 0): "../../../assets/xs/mgxs_2a_Bo6Tf12Tm5.h5",
    (0, 1, 1): "../../../assets/xs/mgxs_2a_Bo6Tf12Tm7.h5",
    (1, 0, 0): "../../../assets/xs/mgxs_2a_Bo13Tf9Tm5.h5",
    (1, 0, 1): "../../../assets/xs/mgxs_2a_Bo13Tf9Tm7.h5",
    (1, 1, 0): "../../../assets/xs/mgxs_2a_Bo13Tf12Tm5.h5",
    (1, 1, 1): "../../../assets/xs/mgxs_2a_Bo13Tf12Tm7.h5",
}

EXPECTED_TOTAL = {
    0: 0.21114941740224946,
    17: 0.29189231171073915,
    79: 0.5685284202920918,
    160: 0.8841928869427658,
    240: 17.371243476724697,
    360: 3.015718391307854,
}
EXPECTED_NU_FISSION = {
    0: 0.14143281494580695,
    17: 0.010381259553111918,
    79: 0.02835739884595183,
    160: 0.13641677142431452,
    240: 0.01201462508018609,
    360: 4.359916837374776,
}
EXPECTED_SCATTER = {
    (0, 0): 0.10446146019372349,
    (1, 0): 0.009966468533361517,
    (13, 10): 0.002056872620449352,
    (17, 10): 0.0011247616068402824,
    (27, 17): 0.004577827673377205,
    (50, 27): 2.179021216595788e-06,
    (51, 50): 0.092119697605361,
    (79, 79): 0.1997720793107284,
    (82, 79): 0.06007157153041619,
    (110, 100): 0.017281619960800534,
    (160, 160): 0.3631568249787166,
    (169, 160): 0.026346106882371392,
    (220, 200): 5.8150365867224216e-05,
    (260, 240): 0.015275496006585076,
    (275, 240): 7.708907915567506e-05,
    (281, 260): 1.2792279664054307e-07,
    (332, 300): 9.907684966585876e-08,
    (360, 360): 0.15125759818503318,
}
EXPECTED_SCATTER_COLUMNS_PER_ROW = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    49, 50, 50, 50, 51, 51, 51, 49, 48, 49, 49, 48, 48, 48, 46, 45,
    46, 45, 43, 42, 38, 37, 38, 38, 38, 37, 36, 36, 37, 36, 37, 38,
    38, 38, 37, 35, 35, 35, 35, 36, 34, 32, 31, 30, 31, 32, 31, 30,
    30, 31, 31, 29, 29, 29, 30, 31, 32, 32, 33, 34, 35, 35, 35, 35,
    34, 30, 26, 26, 25, 26, 26, 27, 27, 27, 28, 26, 26, 27, 28, 28,
    29, 28, 27, 25, 24, 22, 18, 18, 18, 19, 19, 20, 21, 22, 23, 23,
    24, 25, 26, 26, 25, 25, 24, 24, 18, 18, 18, 19, 19, 21, 22, 22,
    21, 24, 24, 25, 23, 24, 24, 24, 24, 23, 23, 23, 21, 20, 20, 20,
    24, 24, 25, 27, 27, 29, 30, 30, 30, 29, 29, 29, 30, 30, 32, 32,
    32, 32, 36, 36, 36, 36, 37, 39, 39, 39, 45, 44, 44, 45, 45, 43,
    41, 33, 41, 37, 37, 35, 34, 34, 34, 34, 35, 35, 34, 32, 36, 36,
    38, 36, 36, 32, 34, 33, 35, 34, 34, 34, 37, 46, 51, 51, 51, 51,
    52, 52, 51, 50, 50, 50, 51, 51, 51, 51, 51, 51, 52, 51, 51, 50,
    53, 53, 54, 53, 54, 54, 54, 54, 54, 55, 55, 57, 57, 55, 57, 57,
    56, 54, 54, 53, 44, 41, 38, 42, 40, 39, 34, 33, 32, 33, 32, 31,
    32, 32, 33, 33, 37, 38, 39, 40, 41, 45, 50, 47, 49, 50, 49, 49,
    50, 51, 51, 51, 52, 50, 50, 49, 49, 50, 50, 49, 48, 48, 49, 49,
    48, 48, 49, 51, 52, 52, 51, 53, 54, 54, 56, 56, 57, 56, 54, 54,
    54, 49, 47, 44, 41, 39, 37, 36, 36, 35, 35, 33, 33, 33, 33, 32,
    32, 30, 30, 30, 28, 28, 27, 26, 25,
]
EXPECTED_PRODUCTION_MATRIX = {
    (0, 0): 4.221782289521071e-06,
    (13, 17): 0.0009623578328127921,
    (17, 79): 0.0006044043852719158,
    (27, 160): 0.002031726376012728,
    (38, 240): 2.7961032773966743e-05,
    (50, 360): 0.0017417129813974105,
}
TOLERANCE = 1.0e-12


def main():
    grid = CartesianGrid(GRID_DATA)
    xs_data = np.empty((2, 2, 2), dtype=object)
    for corner, file_name in CORNER_FILES.items():
        xs = MultiGroupXS()
        xs.LoadFromOpenMC(file_name, "fuel", 294.0)
        xs_data[corner] = xs

    interpolator = LinearInterpolator(
        grid,
        xs_data
    )
    result = interpolator.Evaluate(STATE_POINT)

    checks = [result.num_groups == 361, result.is_fissionable]
    for group, expected in EXPECTED_TOTAL.items():
        actual = float(result.sigma_t[group])
        checks.append(abs(actual - expected) < TOLERANCE)

    for group, expected in EXPECTED_NU_FISSION.items():
        actual = float(result.nu_sigma_f[group])
        checks.append(abs(actual - expected) < TOLERANCE)

    scatter_matrix = result.GetTransferMatrix(0)
    checks.append(scatter_matrix.num_columns_per_row == EXPECTED_SCATTER_COLUMNS_PER_ROW)
    for (row, column), expected in EXPECTED_SCATTER.items():
        actual = scatter_matrix.GetValueIJ(row, column)
        checks.append(abs(actual - expected) < TOLERANCE)

    production_matrix = result.production_matrix
    for (row, column), expected in EXPECTED_PRODUCTION_MATRIX.items():
        actual = production_matrix[row][column]
        checks.append(abs(actual - expected) < TOLERANCE)

    print(f"PASS {int(all(checks))}")


if __name__ == "__main__":
    main()
