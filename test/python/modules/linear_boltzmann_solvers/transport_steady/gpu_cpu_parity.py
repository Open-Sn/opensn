#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""Compare complete CPU and GPU scalar-flux vectors for both GPU sweep types."""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def solve(sweep_type, use_gpus):
    nodes_x = [i / 6.0 for i in range(7)]
    nodes_y = [i / 5.0 for i in range(6)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes_x, nodes_y]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.4)
    source = VolumetricSource(block_ids=[0], group_strength=[1.0])
    quadrature = GLCProductQuadrature2DXY(
        n_polar=2,
        n_azimuthal=4,
        scattering_order=0,
    )

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source],
        boundary_conditions=[
            {"name": "xmin", "type": "isotropic", "group_strength": [0.25]},
            {"name": "ymin", "type": "reflecting"},
        ],
        options={"max_ags_iterations": 1},
        sweep_type=sweep_type,
        use_gpus=use_gpus,
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    return list(problem.GetPhiNewLocal()), solver.ComputeBalanceTable()["balance"]


def compare(reference, candidate):
    if len(reference) != len(candidate):
        raise RuntimeError("CPU and GPU scalar-flux vectors have different sizes.")

    max_abs = max(abs(a - b) for a, b in zip(reference, candidate))
    diff_norm = math.sqrt(sum((a - b) ** 2 for a, b in zip(reference, candidate)))
    ref_norm = math.sqrt(sum(value * value for value in reference))
    return max_abs, diff_norm / max(ref_norm, 1.0e-300)


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Incorrect number of processors. Expected 1 processor but got {size}.")

    for sweep_type in ("AAH", "CBC"):
        cpu_phi, cpu_balance = solve(sweep_type, False)
        gpu_phi, gpu_balance = solve(sweep_type, True)
        max_abs, relative_l2 = compare(cpu_phi, gpu_phi)

        if rank == 0:
            print(f"{sweep_type}_GPU_CPU_MAX_ABS_DIFF={max_abs:.16e}")
            print(f"{sweep_type}_GPU_CPU_REL_L2_DIFF={relative_l2:.16e}")
            print(f"{sweep_type}_GPU_CPU_BALANCE_DIFF={abs(cpu_balance - gpu_balance):.16e}")
