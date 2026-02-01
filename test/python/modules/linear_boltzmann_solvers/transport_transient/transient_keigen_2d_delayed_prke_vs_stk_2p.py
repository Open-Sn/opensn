#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D delayed transient: homogeneous step xs vs PRKE (2 precursors).

1-group, 2 precursors, reflecting boundaries (infinite-medium). A step
to a supercritical xs is applied at t=0. Space-time kinetics should follow
PRKE for this homogeneous case.

PRKE_STK_PASS is 1 if the space-time fission-production ratio matches PRKE within 2%
for t<=0.2.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        TransientSolver,
    )
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator


def read_block_values(path, block_name):
    begin = f"{block_name}_BEGIN"
    end = f"{block_name}_END"
    values = []
    in_block = False
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line == begin:
                in_block = True
                continue
            if line == end:
                in_block = False
                continue
            if in_block:
                parts = line.split()
                if len(parts) >= 2:
                    values.append(float(parts[1]))
    if not values:
        raise RuntimeError(f"Failed to find {block_name} in {path}")
    return values


def solve_linear(A, b):
    n = len(b)
    a = [row[:] for row in A]
    x = b[:]
    for i in range(n):
        pivot = i
        for r in range(i + 1, n):
            if abs(a[r][i]) > abs(a[pivot][i]):
                pivot = r
        if abs(a[pivot][i]) < 1.0e-14:
            raise RuntimeError("Singular system in PRKE solve")
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            x[i], x[pivot] = x[pivot], x[i]
        piv = a[i][i]
        for j in range(i, n):
            a[i][j] /= piv
        x[i] /= piv
        for r in range(n):
            if r == i:
                continue
            factor = a[r][i]
            if factor == 0.0:
                continue
            for j in range(i, n):
                a[r][j] -= factor * a[i][j]
            x[r] -= factor * x[i]
    return x


def prke_step(phi, C, dt, beta, lambdas, rho, Lambda):
    m = len(lambdas)
    beta_total = sum(beta)
    size = 1 + m
    A = [[0.0 for _ in range(size)] for _ in range(size)]
    b = [0.0 for _ in range(size)]

    A[0][0] = 1.0 - dt * (rho - beta_total) / Lambda
    for i in range(m):
        A[0][1 + i] = -dt * lambdas[i]

    b[0] = phi

    for i in range(m):
        A[1 + i][0] = -dt * (beta[i] / Lambda)
        A[1 + i][1 + i] = 1.0 + dt * lambdas[i]
        b[1 + i] = C[i]

    x = solve_linear(A, b)
    return x[0], x[1:]


if __name__ == "__main__":
    dx = 6.0 / 6
    nodes = [i * dx for i in range(6 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_delayed_crit_2p.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_delayed_super_2p.cxs"))

    pquad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)

    num_groups = 1
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    keigen = PowerIterationKEigenSolver(problem=phys, max_iters=200, k_tol=1.0e-10)
    keigen.Initialize()
    keigen.Execute()

    solver = TransientSolver(problem=phys, initial_state="existing")
    solver.Initialize()

    # Apply homogeneous perturbation at t=0
    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])

    sigma_a = xs_super.sigma_a[0]
    nu_sigma_f = xs_super.nu_sigma_f[0]
    inv_vel = xs_super.inv_velocity[0]
    v = 1.0 / inv_vel

    k_eff = nu_sigma_f / sigma_a
    rho = (k_eff - 1.0) / k_eff
    Lambda = 1.0 / (v * nu_sigma_f)

    beta = read_block_values(
        os.path.join(os.path.dirname(__file__), "xs1g_delayed_super_2p.cxs"),
        "PRECURSOR_FRACTIONAL_YIELDS",
    )
    lambdas = read_block_values(
        os.path.join(os.path.dirname(__file__), "xs1g_delayed_super_2p.cxs"),
        "PRECURSOR_DECAY_CONSTANTS",
    )

    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    fp0 = phys.ComputeFissionProduction("new")

    # PRKE initial conditions for steady state
    phi = 1.0
    C = [beta[i] / (Lambda * lambdas[i]) for i in range(len(lambdas))]

    t_end = 0.2
    rel_tol = 2.0e-2
    ok = True
    while phys.GetTime() < t_end:
        solver.Advance()
        fp_new = phys.ComputeFissionProduction("new")
        t_to = phys.GetTime()

        phi, C = prke_step(phi, C, dt, beta, lambdas, rho, Lambda)
        ratio_num = fp_new / fp0
        ratio_prke = phi
        if abs(ratio_num - ratio_prke) > rel_tol * ratio_prke:
            ok = False

    if rank == 0:
        print(f"PRKE_STK_PASS {1 if ok else 0}")
