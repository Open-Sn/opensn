#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D delayed transient k-eigen with semi-analytic 1-precursor kinetics.

Validate delayed-neutron coupling and precursor update against the closed-form
1-precursor point-kinetics solution for a reactivity step.

1-group, 1 precursor. Point kinetics: dphi/dt = ((rho - beta)/Lambda) * phi +
lambda * C dC/dt   = (beta/Lambda) * phi - lambda * C For a step to rho>0 with
phi(0)=1 and C(0)=beta/(Lambda*lambda), the solution is a sum of two
exponentials with eigenvalues w1,w2. The helper delayed_phi_ratio implements
that exact form.

ANALYTIC_PASS is 1 if |FP_ratio - delayed_phi_ratio(t)| < 2% for all steps up to
t=0.2. Parameters: beta=0.0065, lambda=0.08, k=1.2 => rho=(k-1)/k=0.166666...,
nu_total=2.0, sigma_f=0.18, Lambda = 1/(v*nu*sigma_f) = 1/0.36 ≈ 2.7778. These
are the inputs to delayed_phi_ratio. ANALYTIC_PASS validates delayed source and
precursor updates against the semi-analytic solution.
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
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator


def delayed_phi_ratio(t, beta, lam, rho, Lambda):
    # Point-kinetics 1-precursor step solution with phi(0)=1, C(0)=beta/(Lambda*lam)
    a = (rho - beta) / Lambda - lam
    b = math.sqrt(((rho - beta) / Lambda + lam) ** 2 + 4.0 * beta * lam / Lambda)
    w1 = 0.5 * (a + b)
    w2 = 0.5 * (a - b)

    # C = (beta/Lambda)/(w+lam) * phi for each mode
    k1 = (beta / Lambda) / (w1 + lam)
    k2 = (beta / Lambda) / (w2 + lam)

    # Solve for coefficients c1, c2 from phi(0)=1 and C(0)=beta/(Lambda*lam)
    c2 = (beta / (Lambda * lam) - k1) / (k2 - k1)
    c1 = 1.0 - c2

    return c1 * math.exp(w1 * t) + c2 * math.exp(w2 * t)


if __name__ == "__main__":
    dx = 8.0 / 4
    nodes = [i * dx for i in range(4 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_delayed_crit_1p.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_delayed_super_1p.cxs"))

    pquad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

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

    solver = TransientSolver(problem=phys, verbose=False, initial_state="existing")
    solver.Initialize()

    # Swap to supercritical XS at t=0
    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])

    # Semi-analytic parameters
    beta = 0.0065
    lam = 0.08
    k = 1.2
    rho = (k - 1.0) / k
    # Lambda based on infinite-medium definition (prompt gen time)
    nu_prompt = 1.987
    nu_delayed = 0.013
    sigma_f = 0.180000
    v = 1.0
    nu_sigma_f = sigma_f * (nu_prompt + nu_delayed)
    Lambda = 1.0 / (v * nu_sigma_f)

    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    fp0 = phys.ComputeFissionProduction("new")

    t_end = 0.2
    rel_tol = 2.0e-2
    ok = True
    if rank == 0:
        print("step time ratio_numeric ratio_analytic")
    step = 0
    while phys.GetTime() < t_end:
        step += 1
        solver.Advance()
        fp_new = phys.ComputeFissionProduction("new")
        t_to = phys.GetTime()
        ratio_num = fp_new / fp0
        ratio_ana = delayed_phi_ratio(t_to, beta, lam, rho, Lambda)

        if rank == 0:
            print(f"{step:4d} {t_to:10.4e} {ratio_num:12.6e} {ratio_ana:12.6e}")
        if abs(ratio_num - ratio_ana) > rel_tol * ratio_ana:
            ok = False

    if rank == 0:
        print(f"ANALYTIC_PASS {1 if ok else 0}")
