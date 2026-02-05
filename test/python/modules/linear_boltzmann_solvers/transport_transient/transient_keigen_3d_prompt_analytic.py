#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D prompt-only transient k-eigen with analytic exponential response.

1-group, reflecting BCs, no delayed neutrons. After a step to supercritical:
dphi/dt = alpha * phi, with alpha = nu*sigma_f - sigma_a. In discrete time with
theta=1, the update ratio per step is r = (tau + sigma_s + nu*sigma_f)/(tau +
sigma_t), where tau = v^{-1}/dt. For small dt, r ≈ exp(alpha dt), giving
phi(t)/phi(0) = exp(alpha t).

ANALYTIC_PASS is 1 if |FP_ratio - exp(alpha t)| < 0.5% for all steps up to
t=0.1. With sigma_t=1.0, sigma_s=0.7 => sigma_a=0.3, nu=2, sigma_f=0.18,
alpha=0.36-0.3=0.06. Thus exp(alpha*0.1)=exp(0.006)≈1.0060 (used implicitly in
the comparison). ANALYTIC_PASS validates the time term and prompt fission source
handling against the analytic solution.
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

if __name__ == "__main__":
    dx = 8.0 / 4
    nodes = [i * dx for i in range(4 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_prompt_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_prompt_super.cxs"))

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
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    keigen = PowerIterationKEigenSolver(problem=phys, max_iters=200, k_tol=1.0e-10)
    keigen.Initialize()
    keigen.Execute()

    solver = TransientSolver(problem=phys, verbose=False, initial_state="existing")
    solver.Initialize()
    phi_old = phys.GetPhiOldLocal()
    phi_new = phys.GetPhiNewLocal()
    if rank == 0:
        print("phi_old[0]", phi_old[0], "phi_new[0]", phi_new[0])
    # Swap to supercritical XS at t=0
    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])
    fp_new = phys.ComputeFissionProduction("new")
    if rank == 0:
        print("fp_new (new)", fp_new)
    # Analytic alpha for prompt-only
    sigma_t = 1.0
    sigma_s = 0.7
    sigma_a = sigma_t - sigma_s
    nu = 2.0
    sigma_f = 0.180000
    alpha = nu * sigma_f - sigma_a

    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    # Use the converged flux from the k-eigen solve as the initial state
    fp0 = phys.ComputeFissionProduction("new")

    if rank == 0:
        print("inv_velocity", xs_super.inv_velocity)
        print("dt", dt, "theta", 1.0)
    tau = xs_super.inv_velocity[0] / dt
    r_expected = (tau + sigma_s + nu * sigma_f) / (tau + sigma_t)
    if rank == 0:
        print("tau", tau, "r_expected", r_expected, "r_expected^11", r_expected ** 11)

    t_end = 0.1
    rel_tol = 5.0e-3
    ok = True
    if rank == 0:
        print("step time ratio_numeric ratio_analytic")
    step = 0
    while phys.GetTime() < t_end:
        step += 1
        t_from = phys.GetTime()
        solver.Advance()
        fp_step_new = phys.ComputeFissionProduction("new")
        fp_step_old = phys.ComputeFissionProduction("old")
        if rank == 0:
            print(
                "fp_new(after step)",
                fp_step_new,
                "fp_old(after step)",
                fp_step_old,
            )
        phi_old = phys.GetPhiOldLocal()
        phi_new = phys.GetPhiNewLocal()
        if rank == 0:
            print("phi_old[0]", phi_old[0], "phi_new[0]", phi_new[0])
        fp_new = fp_step_new
        t_to = phys.GetTime()
        ratio_num = fp_new / fp0
        ratio_ana = math.exp(alpha * t_to)

        if rank == 0:
            print(f"{step:4d} {t_to:10.4e} {ratio_num:12.6e} {ratio_ana:12.6e}")
        if abs(ratio_num - ratio_ana) > rel_tol * ratio_ana:
            ok = False

    if rank == 0:
        print(f"ANALYTIC_PASS {1 if ok else 0}")
