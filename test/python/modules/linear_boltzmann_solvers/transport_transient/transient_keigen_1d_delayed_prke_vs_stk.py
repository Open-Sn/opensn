#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D delayed transient: homogeneous step xs vs PRKE.

Validate space-time kinetics against point-reactor kinetics (PRKE) for a
homogeneous perturbation in a homogeneous system.

1-group, 1 precursor, reflecting boundaries (infinite-medium). A step
to a supercritical xs is applied at t=0. Space-time kinetics should follow
PRKE for this homogeneous case.

PRKE_STK_PASS is 1 if the space-time fission-production ratio matches PRKE within 2%
for t<=0.2.
"""

import math
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
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator


def read_precursor_value(path, block_name):
    begin = f"{block_name}_BEGIN"
    end = f"{block_name}_END"
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
                    return float(parts[1])
    raise RuntimeError(f"Failed to find {block_name} in {path}")


def prke_phi_ratio(t, beta, lam, rho, Lambda):
    # Point-kinetics 1-precursor step solution with phi(0)=1, C(0)=beta/(Lambda*lam)
    a = (rho - beta) / Lambda - lam
    b = math.sqrt(((rho - beta) / Lambda + lam) ** 2 + 4.0 * beta * lam / Lambda)
    w1 = 0.5 * (a + b)
    w2 = 0.5 * (a - b)

    k1 = (beta / Lambda) / (w1 + lam)
    k2 = (beta / Lambda) / (w2 + lam)

    c2 = (beta / (Lambda * lam) - k1) / (k2 - k1)
    c1 = 1.0 - c2

    return c1 * math.exp(w1 * t) + c2 * math.exp(w2 * t)


if __name__ == "__main__":
    dx = 8.0 / 40
    nodes = [i * dx for i in range(40 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_delayed_crit_1p.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_delayed_super_1p.cxs"))

    pquad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

    num_groups = 1
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
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

    phys.SetTimeDependentMode()

    solver = TransientSolver(problem=phys, initial_state="existing")
    solver.Initialize()

    # Apply homogeneous perturbation at t=0
    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])

    # PRKE parameters from xs
    sigma_a = xs_super.sigma_a[0]
    nu_sigma_f = xs_super.nu_sigma_f[0]
    inv_vel = xs_super.inv_velocity[0]
    v = 1.0 / inv_vel

    # k_inf = nu_sigma_f / sigma_a for 1-group infinite medium
    k_eff = nu_sigma_f / sigma_a
    rho = (k_eff - 1.0) / k_eff

    beta = read_precursor_value(
        os.path.join(os.path.dirname(__file__), "xs1g_delayed_super_1p.cxs"),
        "PRECURSOR_FRACTIONAL_YIELDS",
    )
    lam = read_precursor_value(
        os.path.join(os.path.dirname(__file__), "xs1g_delayed_super_1p.cxs"),
        "PRECURSOR_DECAY_CONSTANTS",
    )
    Lambda = 1.0 / (v * nu_sigma_f)

    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    fp0 = phys.ComputeFissionProduction("new")

    t_end = 0.2
    rel_tol = 2.0e-2
    ok = True
    step = 0
    while phys.GetTime() < t_end:
        step += 1
        solver.Advance()
        fp_new = phys.ComputeFissionProduction("new")
        t_to = phys.GetTime()

        ratio_num = fp_new / fp0
        ratio_prke = prke_phi_ratio(t_to, beta, lam, rho, Lambda)
        if abs(ratio_num - ratio_prke) > rel_tol * ratio_prke:
            ok = False

    if rank == 0:
        print(f"PRKE_STK_PASS {1 if ok else 0}")
