#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.math import Vector3

if __name__ == "__main__":

    # Setup mesh
    nodes = []
    N = 100
    L = 30.0
    xmin = 0.0
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()

    # Cross-section data
    num_groups = 168
    grid.SetUniformBlockID(0)
    xs_3_170 = MultiGroupXS()
    xs_3_170.LoadFromOpenSn("xs_168g.xs")

    # Angular quadrature
    pquad = GLProductQuadrature1DSlab(n_polar=80, scattering_order=5)

    # Create solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            }
        ],
        xs_map=[
            {
                "block_ids": [0],
                "xs": xs_3_170
            }
        ],
        scattering_order=0,
        options={
            "wrong-parameter": True,
            "boundary_conditions": [
            ],
        }
    )
    ss_solver = SteadyStateSolver(problem=phys)
