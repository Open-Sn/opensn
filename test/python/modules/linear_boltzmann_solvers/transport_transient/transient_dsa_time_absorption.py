#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Diffusion synthetic acceleration in transient solves.

The theta scheme adds the time absorption tau_g = 1/(v_g theta dt) to the total
cross section of the transport operator, so the WGDSA and TGDSA diffusion
operators must include it: sigma_r + tau and D = 1/(3 (sigma_tr + tau)).
Without it, for sigma_s/sigma_a >> 1 and tau >~ sigma_a, source iteration with
DSA over-corrects the flat error mode and diverges (flat-mode error factor
c + (sigma_s/sigma_a)(c - 1) with c = sigma_s/(sigma_t + tau)).

Two problems with vacuum boundaries are stepped with backward Euler at small and
changing time steps, which rebuilds the DSA operators:
  1. one group, c = 0.99, Richardson with WGDSA;
  2. two groups with strong upscatter, Richardson with WGDSA and TGDSA, with a
     cross-section swap (rebuilding the two-grid data) partway through.
Each is compared against GMRES without acceleration at a tight tolerance. On a
single rank the WGDSA and TGDSA diffusion solves use direct factorizations, so
case 2 also covers rebuilding two coexisting direct diffusion solvers.

Check:
  DSA_TRANSIENT_MAX_REL_DIFF: max relative difference in the scalar flux over all
  steps and cases. The accelerated runs must converge within the iteration
  limit, so it is at the level of the iteration tolerances.
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return comm.allreduce(value, op=_mpi_ops[op])

    def MPIBarrier():
        comm.Barrier()

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver


def write_xs(path, sigma_t, transfer, inv_velocity):
    # transfer[g_to][g_from]
    if rank == 0:
        with open(path, "w") as f:
            num_groups = len(sigma_t)
            f.write(f"NUM_GROUPS {num_groups}\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
            for g, val in enumerate(sigma_t):
                f.write(f"{g} {val}\n")
            f.write("SIGMA_T_END\n\nINV_VELOCITY_BEGIN\n")
            for g, val in enumerate(inv_velocity):
                f.write(f"{g} {val}\n")
            f.write("INV_VELOCITY_END\n\nTRANSFER_MOMENTS_BEGIN\n")
            for g_to in range(num_groups):
                for g_from in range(num_groups):
                    if transfer[g_to][g_from] != 0.0:
                        f.write(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {transfer[g_to][g_from]}\n")
            f.write("TRANSFER_MOMENTS_END\n")
    MPIBarrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    return xs


def run(grid, xs_list, num_groups, method, wgdsa, tgdsa, steps, swap_step):
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4,
                                                               scattering_order=0),
                "inner_linear_method": method,
                "l_abs_tol": 1.0e-10,
                "l_max_its": 100 if method == "classic_richardson" else 500,
                "apply_wgdsa": wgdsa,
                "apply_tgdsa": tgdsa,
                "wgdsa_l_abs_tol": 1.0e-12,
                "tgdsa_l_abs_tol": 1.0e-12,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs_list[0]}],
        boundary_conditions=[
            {"name": name, "type": "vacuum"} for name in ("xmin", "xmax", "ymin", "ymax")
        ],
        volumetric_sources=[
            VolumetricSource(block_ids=[0], group_strength=[1.0] + [0.0] * (num_groups - 1))
        ],
        options={"save_angular_flux": True, "verbose_inner_iterations": False},
    )
    problem.SetTimeDependentMode()
    solver = TransientSolver(problem=problem, dt=steps[0], stop_time=1.0e9, theta=1.0,
                             initial_state="zero", verbose=False)
    solver.Initialize()
    history = []
    for n, dt in enumerate(steps):
        if n == swap_step:
            problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_list[1]}])
        solver.SetTimeStep(dt)
        solver.Advance()
        history.append(np.array(problem.GetPhiNewLocal(), copy=True))
    return history


if __name__ == "__main__":
    nodes = [i * 10.0 / 8 for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs_1g = write_xs("dsa_td_1g.xs", [1.0], [[0.99]], [1.0])
    xs_2g_a = write_xs("dsa_td_2g_a.xs", [1.0, 2.0], [[0.5, 0.09], [0.45, 1.9]],
                       [1.0, 3.0])
    xs_2g_b = write_xs("dsa_td_2g_b.xs", [1.0, 2.0], [[0.5, 0.19], [0.45, 1.8]],
                       [1.0, 3.0])

    steps = [1.0, 0.1, 0.01, 0.3]
    cases = [
        ([xs_1g, xs_1g], 1, True, False, None),
        ([xs_2g_a, xs_2g_b], 2, True, True, 2),
    ]
    max_diff = 0.0
    for xs_list, num_groups, wgdsa, tgdsa, swap_step in cases:
        ref = run(grid, xs_list, num_groups, "petsc_gmres", False, False, steps, swap_step)
        acc = run(grid, xs_list, num_groups, "classic_richardson", wgdsa, tgdsa, steps, swap_step)
        for phi_ref, phi_acc in zip(ref, acc):
            local_diff = float(np.max(np.abs(phi_acc - phi_ref))) if phi_ref.size else 0.0
            local_scale = float(np.max(np.abs(phi_ref))) if phi_ref.size else 0.0
            diff = MPIAllReduce(local_diff, "max") / MPIAllReduce(local_scale, "max")
            max_diff = max(max_diff, diff if np.isfinite(diff) else 1.0e300)

    MPIBarrier()
    if rank == 0:
        for name in ("dsa_td_1g", "dsa_td_2g_a", "dsa_td_2g_b"):
            os.remove(f"{name}.xs")
        print(f"DSA_TRANSIENT_MAX_REL_DIFF {max_diff:.6e}")
