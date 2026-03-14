#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2D unstructured multi-material equivalence check (density API vs scaled xs).

A two-block problem is solved two ways:
1. Same microscopic xs in both blocks with by-block density scaling.
2. Explicitly scaled macroscopic xs per block with unit density.
Both cases are identical, so scalar flux vectors should match.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver


def make_grid():
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/triangle_mesh2x2_cuts.obj",
        partitioner=KBAGraphPartitioner(nx=2, ny=2, nz=1, xcuts=[0.0], ycuts=[0.0]),
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromFunction(lambda c, old: 1 if c.x > 1.0 else old)
    return grid


def solve_case(use_density_api):
    grid = make_grid()

    quad = GLCProductQuadrature2DXY(n_polar=6, n_azimuthal=24, scattering_order=0)
    src = VolumetricSource(block_ids=[0, 1], group_strength=[1.0])

    if use_density_api:
        xs_micro = MultiGroupXS()
        xs_micro.CreateSimpleOneGroup(1.0, 0.0)
        xs_map = [
            {"block_ids": [0], "xs": xs_micro},
            {"block_ids": [1], "xs": xs_micro},
        ]
        density = {
            "default_density": 1.0,
            "by_block": [{"block_ids": [1], "density": 2.5}],
        }
    else:
        xs0 = MultiGroupXS()
        xs0.CreateSimpleOneGroup(1.0, 0.0)
        xs1 = MultiGroupXS()
        xs1.CreateSimpleOneGroup(2.5, 0.0)
        xs_map = [
            {"block_ids": [0], "xs": xs0},
            {"block_ids": [1], "xs": xs1},
        ]
        density = {"default_density": 1.0}

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
        xs_map=xs_map,
        volumetric_sources=[src],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
        ],
        density=density,
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    return list(problem.GetPhiNewLocal())


if __name__ == "__main__":
    phi_api = solve_case(True)
    phi_ref = solve_case(False)

    diff2 = sum((a - b) * (a - b) for a, b in zip(phi_api, phi_ref))
    ref2 = sum(b * b for b in phi_ref)
    rel = (diff2 / ref2) ** 0.5 if ref2 > 0.0 else 0.0
    pass_flag = 1 if rel < 1.0e-10 else 0

    if rank == 0:
        print(f"DENSITY_MM_BAL_REL={rel:.8e}")
        print(f"DENSITY_MM_BAL_PASS {pass_flag}")
