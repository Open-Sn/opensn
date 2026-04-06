#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D delayed transient with XS swaps across fissionability and precursor counts.

This regression exercises runtime XS-map swaps that:
1. change from a delayed fissionable material with 2 precursor families
2. to a non-fissionable pure absorber with 0 precursor families
3. to a delayed fissionable material with 1 precursor family

The geometry is homogeneous with reflecting boundaries, so the transport solution
reduces to an infinite-medium kinetics problem. We compare the scalar-flux ratio
against an exact Backward Euler update for the corresponding reduced model:

- During the non-fissionable interval, phi_{n+1} = phi_n / (1 + v sigma_t dt)
- During the 1-precursor interval, [phi, C]^T is advanced by the exact 2x2
  Backward Euler kinetics solve with C initialized to zero after the 0-precursor swap

PRECURSOR_SWAP_PASS is 1 if the numeric scalar-flux ratio matches the piecewise
reference within 5.0e-6 relative error at every step.
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


def backward_euler_prompt_decay(phi_n, dt, sigma_t, velocity):
    return phi_n / (1.0 + velocity * sigma_t * dt)


def backward_euler_one_precursor(phi_n, c_n, dt, reactivity, beta, decay_const, generation_time):
    a_coeff = (reactivity - beta) / generation_time
    b_coeff = beta / generation_time

    A = 1.0 - dt * a_coeff
    D = 1.0 + dt * decay_const
    det = A * D - (dt * decay_const) * (dt * b_coeff)

    phi_np1 = (D * phi_n + dt * decay_const * c_n) / det
    c_np1 = (dt * b_coeff * phi_n + A * c_n) / det
    return phi_np1, c_np1


if __name__ == "__main__":
    n_cells = 40
    length = 8.0
    dx = length / n_cells
    nodes = [i * dx for i in range(n_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    base_dir = os.path.dirname(__file__)
    xs_crit_2p_path = os.path.join(base_dir, "xs1g_delayed_crit_2p.cxs")
    xs_super_1p_path = os.path.join(base_dir, "xs1g_delayed_super_1p.cxs")

    xs_crit_2p = MultiGroupXS()
    xs_crit_2p.LoadFromOpenSn(xs_crit_2p_path)

    xs_absorber = MultiGroupXS()
    xs_absorber.CreateSimpleOneGroup(1.5, 0.0, 1.0)

    xs_super_1p = MultiGroupXS()
    xs_super_1p.LoadFromOpenSn(xs_super_1p_path)

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
        xs_map=[{"block_ids": [0], "xs": xs_crit_2p}],
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

    keigen = PowerIterationKEigenSolver(problem=problem, max_iters=200, k_tol=1.0e-10)
    keigen.Initialize()
    keigen.Execute()

    problem.SetTimeDependentMode()

    solver = TransientSolver(problem=problem, verbose=False, initial_state="existing")
    solver.Initialize()

    phi_initial = max_phi(problem)
    if phi_initial <= 0.0:
        raise RuntimeError("Initial scalar flux must be positive.")

    dt_absorb = 5.0e-2
    dt_fission = 1.0e-2
    n_absorb_steps = 3
    n_fission_steps = 10
    rel_tol = 5.0e-6

    ok = True
    phi_ref = 1.0
    c_ref = 0.0

    absorber_sigma_t = xs_absorber.sigma_t[0]
    absorber_velocity = 1.0 / xs_absorber.inv_velocity[0]

    beta = read_precursor_values(xs_super_1p_path, "PRECURSOR_FRACTIONAL_YIELDS")[0]
    decay_const = read_precursor_values(xs_super_1p_path, "PRECURSOR_DECAY_CONSTANTS")[0]
    sigma_a = xs_super_1p.sigma_a[0]
    nu_sigma_f = xs_super_1p.nu_sigma_f[0]
    velocity = 1.0 / xs_super_1p.inv_velocity[0]
    k_eff = nu_sigma_f / sigma_a
    reactivity = (k_eff - 1.0) / k_eff
    generation_time = 1.0 / (velocity * nu_sigma_f)

    if rank == 0:
        print("phase step time ratio_numeric ratio_reference")

    problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_absorber}])
    solver.SetTimeStep(dt_absorb)
    solver.SetTheta(1.0)

    for step in range(1, n_absorb_steps + 1):
        solver.Advance()
        phi_ref = backward_euler_prompt_decay(
            phi_ref, dt_absorb, absorber_sigma_t, absorber_velocity
        )
        phi_num = max_phi(problem)
        ratio_num = phi_num / phi_initial
        time_now = problem.GetTime()
        if rank == 0:
            print(f"absorb {step:4d} {time_now:10.4e} {ratio_num:12.6e} {phi_ref:12.6e}")
        if abs(ratio_num - phi_ref) > rel_tol * max(phi_ref, 1.0e-12):
            ok = False

    problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super_1p}])
    solver.SetTimeStep(dt_fission)
    c_ref = 0.0

    for step in range(1, n_fission_steps + 1):
        solver.Advance()
        phi_ref, c_ref = backward_euler_one_precursor(
            phi_ref, c_ref, dt_fission, reactivity, beta, decay_const, generation_time
        )
        phi_num = max_phi(problem)
        ratio_num = phi_num / phi_initial
        time_now = problem.GetTime()
        if rank == 0:
            print(f"fission {step:3d} {time_now:10.4e} {ratio_num:12.6e} {phi_ref:12.6e}")
        if abs(ratio_num - phi_ref) > rel_tol * max(phi_ref, 1.0e-12):
            ok = False

    if rank == 0:
        print(f"PRECURSOR_SWAP_PASS {1 if ok else 0}")
