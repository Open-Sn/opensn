#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D LinearBSolver test of a block of graphite with an air cavity. DSA and TG
SDM: PWLD
Test: WGS groups [0-62] Iteration    28 Residual 6.74299e-07 CONVERGED
and   WGS groups [63-167] Iteration    39 Residual 8.73816e-07 CONVERGED
"""

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
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    num_procs = 4

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    nodes = []
    N = 1000
    L = 100
    xmin = -L / 2
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    vol1 = RPPLogicalVolume(infx=True, infy=True, zmin=-10.0, zmax=10.0)
    grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    num_groups = 168
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("xs_graphite_pure.xs")
    xs_air = MultiGroupXS()
    xs_air.LoadFromOpenSn("xs_air50RH.xs")

    strength = [0.0 for _ in range(num_groups)]
    strength[0] = 1.0
    mg_src = VolumetricSource(block_ids=[0], group_strength=strength)

    # Setup Physics
    pquad = GLProductQuadrature1DSlab(n_polar=4)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, 62),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 1000,
                "gmres_restart_interval": 30,
                "apply_wgdsa": True,
                "wgdsa_l_abs_tol": 1.0e-2,
            },
            {
                "groups_from_to": (63, num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 1000,
                "gmres_restart_interval": 30,
                "apply_wgdsa": True,
                "apply_tgdsa": True,
                "wgdsa_l_abs_tol": 1.0e-2,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_graphite},
            {"block_ids": [1], "xs": xs_air},
        ],
        options={"scattering_order": 1,
                 "max_ags_iterations": 1,
                 "volumetric_sources": [mg_src]
                 },
    )

    # Initialize and Execute Solver
    ss_solver = SteadyStateSolver(lbs_problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()
