#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D fission transient with an XS swap from 0 to 2 precursor families.

This regression starts from a homogeneous prompt-only fissionable critical
k-eigen state with zero precursor families, then swaps to a delayed
2-precursor fissionable material before advancing in time.

The geometry is homogeneous with reflecting boundaries, so after the swap the
space-time transport solution reduces to homogeneous point kinetics with two
precursor families. The reference model uses the exact Backward Euler update
for the same kinetics equations and initial condition C_j(0)=0.

PRECURSOR_INTRO_PASS is 1 if the scalar-flux ratio matches the piecewise
Backward Euler reference within 5.0e-6 relative error for all tested time steps.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        TransientSolver,
    )
    from pyopensn.xs import MultiGroupXS


def read_precursor_values(path, block_name):
    begin = f"{block_name}_BEGIN"
    end = f"{block_name}_END"
    in_block = False
    values = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line == begin:
                in_block = True
                continue
            if line == end:
                break
            if in_block:
                parts = line.split()
                if len(parts) >= 2:
                    values.append(float(parts[1]))
    if not values:
        raise RuntimeError(f"Failed to find values for {block_name} in {path}")
    return values


def max_phi(problem):
    fflist = problem.GetScalarFluxFieldFunction()
    monitor_volume = RPPLogicalVolume(infx=True, infy=True, infz=True)
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType("max")
    field_interp.SetLogicalVolume(monitor_volume)
    field_interp.AddFieldFunction(fflist[0])
    field_interp.Execute()
    return field_interp.GetValue()


def backward_euler_two_precursors(phi_n, c_n, dt, reactivity, betas, decay_consts, generation_time):
    beta_total = sum(betas)
    a_coeff = (reactivity - beta_total) / generation_time

    denominator = 1.0 - dt * a_coeff
    numerator = phi_n
    new_precursors = []

    for beta_j, lam_j, c_j_n in zip(betas, decay_consts, c_n):
        source_coeff = beta_j / generation_time
        denom_j = 1.0 + dt * lam_j
        denominator -= (dt * dt * lam_j * source_coeff) / denom_j
        numerator += (dt * lam_j * c_j_n) / denom_j

    phi_np1 = numerator / denominator

    for beta_j, lam_j, c_j_n in zip(betas, decay_consts, c_n):
        source_coeff = beta_j / generation_time
        c_np1 = (c_j_n + dt * source_coeff * phi_np1) / (1.0 + dt * lam_j)
        new_precursors.append(c_np1)

    return phi_np1, new_precursors


if __name__ == "__main__":
    n_cells = 40
    length = 8.0
    dx = length / n_cells
    nodes = [i * dx for i in range(n_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    base_dir = os.path.dirname(__file__)
    xs_prompt_crit_path = os.path.join(base_dir, "xs1g_prompt_crit.cxs")
    xs_super_2p_path = os.path.join(base_dir, "xs1g_delayed_super_2p.cxs")

    xs_prompt_crit = MultiGroupXS()
    xs_prompt_crit.LoadFromOpenSn(xs_prompt_crit_path)

    xs_super_2p = MultiGroupXS()
    xs_super_2p.LoadFromOpenSn(xs_super_2p_path)

    pquad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs_prompt_crit}],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )

    keigen = PowerIterationKEigenSolver(problem=problem, max_iters=200, k_tol=1.0e-10)
    keigen.Initialize()
    keigen.Execute()

    phi_initial = max_phi(problem)
    if phi_initial <= 0.0:
        raise RuntimeError("Initial scalar flux must be positive.")

    problem.SetTimeDependentMode()

    solver = TransientSolver(problem=problem, verbose=False, initial_state="existing")
    solver.Initialize()

    problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super_2p}])

    betas = read_precursor_values(xs_super_2p_path, "PRECURSOR_FRACTIONAL_YIELDS")
    decay_consts = read_precursor_values(xs_super_2p_path, "PRECURSOR_DECAY_CONSTANTS")
    sigma_a = xs_super_2p.sigma_a[0]
    nu_sigma_f = xs_super_2p.nu_sigma_f[0]
    velocity = 1.0 / xs_super_2p.inv_velocity[0]
    k_eff = nu_sigma_f / sigma_a
    reactivity = (k_eff - 1.0) / k_eff
    generation_time = 1.0 / (velocity * nu_sigma_f)

    dt = 1.0e-2
    n_steps = 20
    rel_tol = 5.0e-6

    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    phi_ref = 1.0
    c_ref = [0.0 for _ in betas]
    ok = True

    if rank == 0:
        print("step time ratio_numeric ratio_reference")

    for step in range(1, n_steps + 1):
        solver.Advance()
        phi_ref, c_ref = backward_euler_two_precursors(
            phi_ref, c_ref, dt, reactivity, betas, decay_consts, generation_time
        )
        phi_num = max_phi(problem)
        ratio_num = phi_num / phi_initial
        time_now = problem.GetTime()
        if rank == 0:
            print(f"{step:4d} {time_now:10.4e} {ratio_num:12.6e} {phi_ref:12.6e}")
        if abs(ratio_num - phi_ref) > rel_tol * max(phi_ref, 1.0e-12):
            ok = False

    if rank == 0:
        print(f"PRECURSOR_INTRO_PASS {1 if ok else 0}")
