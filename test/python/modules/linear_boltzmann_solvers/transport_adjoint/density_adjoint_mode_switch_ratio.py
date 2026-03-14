#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Forward-adjoint-forward mode-switch consistency with density scaling.

This test checks two things:
1. Repeating forward solve after an adjoint solve reproduces the original forward flux.
2. For a reflected one-group absorber/source system, mean forward flux scales as
   phi proportional to 1/rho, so doubling density should halve the forward mean.
Expected solution: The repeat difference ~= 0 and mean(rho=2)/mean(rho=1) ~= 0.5.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.Get_rank()
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver


def run_density_case(rho):
    nodes = [i / 24.0 for i in range(25)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)

    forward_source = VolumetricSource(block_ids=[0], group_strength=[1.0])
    adjoint_source = VolumetricSource(block_ids=[0], group_strength=[0.7])

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": GLProductQuadrature1DSlab(n_polar=6, scattering_order=0),
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 300,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[forward_source],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        density={"default_density": rho},
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()

    solver.Execute()
    phi_fwd_1 = list(problem.GetPhiNewLocal())

    problem.SetAdjoint(True)
    problem.SetBoundaryOptions(clear_boundary_conditions=True)
    problem.SetBoundaryOptions(
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ]
    )
    problem.SetVolumetricSources(volumetric_sources=[adjoint_source])
    solver.Execute()

    problem.SetAdjoint(False)
    problem.SetBoundaryOptions(clear_boundary_conditions=True)
    problem.SetBoundaryOptions(
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ]
    )
    problem.SetVolumetricSources(volumetric_sources=[forward_source])
    solver.Execute()
    phi_fwd_2 = list(problem.GetPhiNewLocal())

    max_abs_diff = max(abs(a - b) for a, b in zip(phi_fwd_1, phi_fwd_2))
    fwd_mean = sum(phi_fwd_2) / len(phi_fwd_2)
    return fwd_mean, max_abs_diff


if __name__ == "__main__":
    mean1, repeat_diff1 = run_density_case(1.0)
    mean2, repeat_diff2 = run_density_case(2.0)

    ratio = mean2 / mean1 if mean1 != 0.0 else 0.0
    repeat_ok = (repeat_diff1 < 1.0e-9) and (repeat_diff2 < 1.0e-9)
    ratio_ok = abs(ratio - 0.5) < 1.0e-2
    pass_flag = 1 if (repeat_ok and ratio_ok) else 0

    if rank == 0:
        print(f"DENSITY_ADJOINT_REPEAT_DIFF_R1={repeat_diff1:.8e}")
        print(f"DENSITY_ADJOINT_REPEAT_DIFF_R2={repeat_diff2:.8e}")
        print(f"DENSITY_ADJOINT_RATIO={ratio:.8e}")
        print(f"DENSITY_ADJOINT_PASS {pass_flag}")
