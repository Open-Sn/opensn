#!/usr/bin/env python3

import numpy as np

if "opensn_console" not in globals():
    from pyopensn.xs import CartesianGrid, LinearInterpolator, MultiGroupXS, XSType


def one_group_total(x, y, z):
    return 2.0 + 0.2 * x - 0.1 * y + 0.05 * z


def one_group_scatter(x, y, z):
    return 0.4 + 0.03 * x + 0.02 * y - 0.01 * z


def make_simple_one_group_xs(x, y, z):
    sigma_t = one_group_total(x, y, z)
    sigma_s = one_group_scatter(x, y, z)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=sigma_s / sigma_t)
    return xs


def main():
    grid_data = [[0.0, 1.5], [-1.0, 2.0], [10.0, 13.0]]
    grid = CartesianGrid(grid_data)

    xs_data = np.empty((2, 2, 2), dtype=object)
    for i, x in enumerate(grid_data[0]):
        for j, y in enumerate(grid_data[1]):
            for k, z in enumerate(grid_data[2]):
                xs_data[i, j, k] = make_simple_one_group_xs(x, y, z)

    interpolator = LinearInterpolator(
        grid,
        xs_data,
        XSType.Total | XSType.Absorption | XSType.Transfer,
    )

    state_point = np.array([0.63, 0.4, 11.7], dtype=float)
    result = interpolator.Evaluate(state_point)

    expected_total = one_group_total(*state_point)
    expected_scatter = one_group_scatter(*state_point)
    expected_absorption = expected_total - expected_scatter

    checks = [
        result.num_groups == 1,
        len(result.sigma_t) == 1,
        len(result.sigma_a) == 1,
        abs(result.sigma_t[0] - expected_total) < 1.0e-12,
        abs(result.sigma_a[0] - expected_absorption) < 1.0e-12,
        abs(result.GetTransferMatrix(0).GetValueIJ(0, 0) - expected_scatter) < 1.0e-12,
    ]

    print(f"interp sigma_t {result.sigma_t[0]:.16e}")
    print(f"interp sigma_a {result.sigma_a[0]:.16e}")
    print(f"interp scatter {result.GetTransferMatrix(0).GetValueIJ(0, 0):.16e}")
    print(f"PASS {int(all(checks))}")


if __name__ == "__main__":
    main()
