#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cross-section scaling test: xs.Combine vs density API.

Case A: microscopic XS scaled using XS.Combine([(xs, scale)]) with density=1.
Case B: unscaled microscopic XS with density=scale.
Both cases represent the same macroscopic model, so scalar flux vectors should match
(relative difference ~= 0).
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver


def solve_case(use_density_api):
    nodes = [i / 12.0 for i in range(13)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs_micro = MultiGroupXS()
    xs_micro.CreateSimpleOneGroup(1.0, 0.0)

    scale = 1.7
    if use_density_api:
        xs = xs_micro
        density = {"default_density": scale}
    else:
        xs = MultiGroupXS()
        xs.Combine([(xs_micro, scale)])
        density = {"default_density": 1.0}

    quad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)
    src = VolumetricSource(block_ids=[0], group_strength=[1.0])

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": quad,
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 200,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=[
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"},
        ],
        density=density,
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return list(problem.GetPhiNewLocal())


if __name__ == "__main__":
    phi_density = solve_case(True)
    phi_combine = solve_case(False)

    diff2 = sum((a - b) * (a - b) for a, b in zip(phi_density, phi_combine))
    ref2 = sum(b * b for b in phi_combine)
    rel = (diff2 / ref2) ** 0.5 if ref2 > 0.0 else 0.0
    pass_flag = 1 if rel < 1.0e-10 else 0

    if rank == 0:
        print(f"DENSITY_COMBINE_REL_L2={rel:.8e}")
        print(f"DENSITY_COMBINE_PASS {pass_flag}")
