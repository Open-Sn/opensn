#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Infinite 172-group problem with balance checks

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
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver

if __name__ == "__main__":

    # Create mesh
    N = 2
    L = 10
    xmin = -L / 2
    dx = L / N
    nodes = [xmin + k * dx for k in range(N + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()

    # Set uniform block IDs
    grid.SetUniformBlockID(0)

    # Set number of energy groups
    num_groups = 172

    # Load cross sections from OpenMC MGXS file
    xs_hdpe = MultiGroupXS()
    xs_hdpe.LoadFromOpenMC("HDPE.h5", "set1", 294.0)

    # Create sources
    strength = [0.0 for _ in range(num_groups)]
    strength[0] = 1.0
    mg_src = VolumetricSource(block_ids=[0], group_strength=strength)

    # Angular quadrature
    pquad = GLCProductQuadrature3DXYZ(4, 8)

    # Solver
    phys = DiscreteOrdinatesSolver(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-9,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_hdpe},
        ],
        options={
            "scattering_order": 0,
            "spatial_discretization": "pwld",
            "save_angular_flux": True,
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
        sweep_type="CBC"
    )
    ss_solver = SteadyStateSolver(lbs_problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Compute particle balance
    phys.ComputeBalance()
