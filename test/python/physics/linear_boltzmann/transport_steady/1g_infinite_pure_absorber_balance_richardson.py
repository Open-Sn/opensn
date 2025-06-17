#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Infinite, 1-group, pure absorber with balance

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver

if __name__ == "__main__":

    # Create Mesh
    nodes = []
    N = 2
    L = 10
    xmin = -L / 2
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    # Add cross sections to materials
    xs1g = MultiGroupXS()
    xs1g.CreateSimpleOneGroup(1.0, 0.0)

    mg_src = VolumetricSource(block_ids=[0], group_strength=[1.0])

    # Angular Quadrature
    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8)

    # LBS block option
    num_groups = 1
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-9,
                "l_max_its": 300,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs1g},
        ],
        options={
            "boundary_conditions": [
                {"name": "xmin", "type": "reflecting"},
                {"name": "xmax", "type": "reflecting"},
                {"name": "ymin", "type": "reflecting"},
                {"name": "ymax", "type": "reflecting"},
                {"name": "zmin", "type": "reflecting"},
                {"name": "zmax", "type": "reflecting"},
            ],
            "volumetric_sources": [mg_src],
        },
    )

    # Initialize and execute solver
    ss_solver = SteadyStateSolver(lbs_problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # compute particle balance
    phys.ComputeBalance()
