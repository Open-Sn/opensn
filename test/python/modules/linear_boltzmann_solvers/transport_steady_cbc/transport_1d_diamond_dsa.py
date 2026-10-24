#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""
1D diamond-difference CBC transport with consistent multigroup DSA
The scalar fluxes are compared with the solution of the
diamond-difference S_N equations, a pure absorber,
and a constant manufactured multigroup solution.
Test: Reference-error, Constant-error, and Balance-error below 1e-9
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.math import Vector3


def reference_flux(nodes, angles, total, self_scatter, downscatter, source, left, right):
    mu, weights = np.polynomial.legendre.leggauss(angles)
    weights /= 2.0
    edges = len(nodes)
    result = []
    for g, sigma_t in enumerate(total):
        matrix = np.zeros((angles * edges, angles * edges))
        rhs = np.zeros(angles * edges)
        for m, cosine in enumerate(mu):
            boundary = 0 if cosine > 0.0 else edges - 1
            row = m * edges
            matrix[row, m * edges + boundary] = 1.0
            rhs[row] = left[g] if cosine > 0.0 else right[g]
            for i, h in enumerate(np.diff(nodes)):
                row = m * edges + i + 1
                matrix[row, m * edges + i] = -cosine + sigma_t * h / 2.0
                matrix[row, m * edges + i + 1] = cosine + sigma_t * h / 2.0
                for n, weight in enumerate(weights):
                    matrix[row, n * edges + i] -= self_scatter[g] * weight * h / 2.0
                    matrix[row, n * edges + i + 1] -= self_scatter[g] * weight * h / 2.0
                rhs[row] = source[g] * h
                if g > 0:
                    rhs[row] += (downscatter[g - 1] * h
                                 * (result[g - 1][i] + result[g - 1][i + 1]) / 2.0)
        angular = np.linalg.solve(matrix, rhs).reshape(angles, edges)
        result.append(weights @ angular)
    return np.asarray(result)


def solve(nodes, angles, xs, source, left, right, accelerated, max_iterations=2000):
    grid = OrthogonalMeshGenerator(node_sets=[list(nodes)]).Execute()
    grid.SetUniformBlockID(0)
    groups = len(source)
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=groups,
        groupsets=[{
            "groups_from_to": (0, groups - 1),
            "angular_quadrature": GLProductQuadrature1DSlab(n_polar=angles, scattering_order=0),
            "inner_linear_method": "classic_richardson",
            "l_abs_tol": 1.0e-11,
            "l_max_its": max_iterations,
            "allow_cycles": False,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=list(source))],
        boundary_conditions=[
            {"name": "zmin", "type": "isotropic", "group_strength": list(left)},
            {"name": "zmax", "type": "isotropic", "group_strength": list(right)},
        ],
        sweep_type="CBC",
        transport_discretization="diamond_difference",
        apply_1d_dsa=accelerated,
        options={"max_ags_iterations": 1},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    if accelerated and max_iterations > 1:
        solver.Execute()
    fields = problem.GetScalarFluxFieldFunction(only_scalar_flux=False)
    values = np.zeros((groups, len(nodes) - 1))
    for g in range(groups):
        for i, z in enumerate((np.asarray(nodes[:-1]) + nodes[1:]) / 2.0):
            point = FieldFunctionInterpolationPoint()
            point.SetPointOfInterest(Vector3(0.0, 0.0, float(z)))
            point.AddFieldFunction(fields[g][0])
            point.Execute()
            values[g, i] = point.GetPointValue()
    return values, abs(solver.ComputeBalanceTable()["balance"])


if __name__ == "__main__":

    num_procs = 3
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    nodes = np.linspace(0.0, 1.0, 13)
    absorber = MultiGroupXS()
    absorber.CreateSimpleOneGroup(1.0, 0.0)
    factor = (1.0 / np.sqrt(3.0) - 1.0 / 24.0) / (1.0 / np.sqrt(3.0) + 1.0 / 24.0)
    edge_flux = 0.5 * factor ** np.arange(len(nodes))
    expected = (edge_flux[:-1] + edge_flux[1:]) / 2.0
    error = 0.0
    balance_error = 0.0
    for accelerated in (False, True):
        values, balance = solve(nodes, 2, absorber, [0.0], [1.0], [0.0], accelerated)
        error = max(error, float(np.max(np.abs(values[0] - expected))))
        balance_error = max(balance_error, balance)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("diamond_downscatter.xs")
    total = [2.25, 2.25, 1.25]
    self_scatter = [1.0, 1.0, 1.0]
    downscatter = [1.0, 1.0]

    constant = np.array([1.0, 2.0, 8.0])
    source = np.array([1.25, 1.5, 0.0])
    values, balance = solve(nodes, 8, xs, source, constant, constant, True)
    constant_error = float(np.max(np.abs(values - constant[:, None])))
    balance_error = max(balance_error, balance)

    nodes = np.array([0.0, 0.03, 0.08, 0.16, 0.30, 0.50, 0.75, 1.0, 1.3, 1.7, 2.2, 3.0, 4.0])
    source = [0.0, 0.0, 0.0]
    left = [0.0, 0.0, 0.0]
    right = [1.0, 0.0, 0.0]
    for angles in (2, 8):
        edges = reference_flux(nodes, angles, total, self_scatter, downscatter,
                               source, left, right)
        expected = (edges[:, :-1] + edges[:, 1:]) / 2.0
        if angles == 2:
            for scale in (1.0, 1.0e-15):
                values, _ = solve(nodes, angles, xs, source, left, [scale * r for r in right],
                                  True, max_iterations=1)
                error = max(error, float(np.max(np.abs(values / scale - expected))))
        for accelerated in (False, True):
            values, balance = solve(nodes, angles, xs, source, left, right, accelerated)
            error = max(error, float(np.max(np.abs(values - expected))))
            balance_error = max(balance_error, balance)

    if rank == 0:
        print(f"Reference-error={error:.12e}")
        print(f"Constant-error={constant_error:.12e}")
        print(f"Balance-error={balance_error:.12e}")
