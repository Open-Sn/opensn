#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 5G infinite-medium k-eigenvalue test with five groupsets and upscatter.

The homogeneous fully reflecting problem has a spatially constant solution.
The reference k is the dominant eigenvalue of inv(T - S) F for the XS data in
simple_5g_upscatter_fissile.xs.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        NonLinearKEigenSolver,
        PowerIterationKEigenSolver,
    )


SIGMA_T = [1.0, 1.0, 1.0, 1.0, 1.0]
SCATTER = [
    [0.20, 0.02, 0.03, 0.01, 0.02],
    [0.04, 0.25, 0.02, 0.03, 0.01],
    [0.06, 0.05, 0.22, 0.04, 0.02],
    [0.02, 0.06, 0.07, 0.24, 0.03],
    [0.01, 0.03, 0.04, 0.05, 0.21],
]
PRODUCTION = [
    [0.275, 0.330, 0.110, 0.055, 0.044],
    [0.125, 0.150, 0.050, 0.025, 0.020],
    [0.060, 0.072, 0.024, 0.012, 0.010],
    [0.030, 0.036, 0.012, 0.006, 0.005],
    [0.010, 0.012, 0.004, 0.002, 0.002],
]
EXPECTED_K = 0.647082071109


def solve_dense(matrix, rhs):
    """Small dense Gaussian solve with partial pivoting."""
    matrix = [row[:] for row in matrix]
    rhs = rhs[:]
    n = len(rhs)
    for i in range(n):
        pivot = max(range(i, n), key=lambda r: abs(matrix[r][i]))
        matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
        rhs[i], rhs[pivot] = rhs[pivot], rhs[i]
        scale = matrix[i][i]
        if abs(scale) < 1.0e-14:
            raise RuntimeError("Singular reference matrix.")
        for j in range(i, n):
            matrix[i][j] /= scale
        rhs[i] /= scale
        for r in range(n):
            if r == i:
                continue
            factor = matrix[r][i]
            for j in range(i, n):
                matrix[r][j] -= factor * matrix[i][j]
            rhs[r] -= factor * rhs[i]
    return rhs


def apply_reference_operator(phi):
    loss = [
        [(SIGMA_T[g] if g == gp else 0.0) - SCATTER[gp][g] for gp in range(5)]
        for g in range(5)
    ]
    fission_source = [
        sum(PRODUCTION[g][gp] * phi[gp] for gp in range(5))
        for g in range(5)
    ]
    return solve_dense(loss, fission_source)


def compute_reference_k():
    phi = [1.0, 1.0, 1.0, 1.0, 1.0]
    k = 1.0
    for _ in range(200):
        new_phi = apply_reference_operator(phi)
        k = sum(new_phi) / sum(phi)
        norm = sum(new_phi)
        phi = [x / norm for x in new_phi]
    return k


def make_problem():
    nodes = [0.0, 1.0, 2.0]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("simple_5g_upscatter_fissile.xs")

    quad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)
    groupsets = []
    for group in range(5):
        groupsets.append(
            {
                "groups_from_to": [group, group],
                "angular_quadrature": quad,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 200,
                "gmres_restart_interval": 30,
            }
        )

    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=5,
        groupsets=groupsets,
        xs_map=[
            {"block_ids": [0], "xs": xs},
        ],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": True,
            "max_ags_iterations": 200,
            "ags_tolerance": 1.0e-10,
        },
    )


def run_solver(name):
    problem = make_problem()
    if name == "PowerIteration":
        solver = PowerIterationKEigenSolver(problem=problem, max_iters=200, k_tol=1.0e-9)
    elif name == "NonLinearKEigen":
        solver = NonLinearKEigenSolver(
            problem=problem,
            nl_max_its=100,
            nl_abs_tol=1.0e-10,
            l_abs_tol=1.0e-10,
            l_max_its=200,
            num_initial_power_iterations=10,
        )
    else:
        raise ValueError(f"Unknown solver {name}.")

    solver.Initialize()
    solver.Execute()
    k = solver.GetEigenvalue()
    error = abs(k - EXPECTED_K)
    if rank == 0:
        print(f"{name} k-eigenvalue: {k:.12f}")
        print(f"{name} k-error: {error:.6e}")
    return name, k, error


if __name__ == "__main__":
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    reference_k = compute_reference_k()
    if abs(reference_k - EXPECTED_K) > 1.0e-12:
        raise RuntimeError(
            f"Internal reference mismatch: computed {reference_k:.16e}, expected {EXPECTED_K:.16e}."
        )

    if rank == 0:
        print(f"Expected k-eigenvalue: {EXPECTED_K:.12f}")

    failures = []
    for solver_name in ("PowerIteration", "NonLinearKEigen"):
        name, k, error = run_solver(solver_name)
        if error > 5.0e-6:
            failures.append(f"{name} k={k:.12f} differs from expected {EXPECTED_K:.12f}.")

    if failures:
        raise RuntimeError(" ; ".join(failures))
