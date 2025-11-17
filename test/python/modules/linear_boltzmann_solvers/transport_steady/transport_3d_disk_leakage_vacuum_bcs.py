#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simple leakage test with semi-analytic solution

Consider a disk with r = 0.2 m, h = 0.1 m, V = pi * r^2 * h = 0.0126 m^3, and
Q = 122.58 particles/s, with all vacuum boundaries.

For a single energy group, pure absorber with sigma_t = 1.0 m^-1, the total leakage can be
approximated as: L = Q*<e^(-sigma_t*s)> where s = mean distance to boundary ~= 0.09m.

Hence, total leakage ~= 122.58*e^(-1.0*0.09) = 112.03 particles/s.

OpenSn result:
Top leakage=36.10
Bottom leakage=36.10
Side leakage=39.75
Total leakage=111.96
"""

import os
import numpy as np
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

if __name__ == "__main__":

    # Mesh
    meshgen = FromFileMeshGenerator(filename="../../../../assets/mesh/disk.msh")
    grid = meshgen.Execute()

    # Materials
    mat1 = MultiGroupXS()
    mat1.CreateSimpleOneGroup(sigma_t=1.0, c=0.0)
    xsecs = [{"block_ids": [1], "xs": mat1}]

    # Volumetric source
    vol_src = VolumetricSource(block_ids=[1], group_strength=[9754.5])

    # Solver
    pquad = GLCProductQuadrature3DXYZ(n_polar=16, n_azimuthal=128, scattering_order=0)

    num_groups = 1
    groupsets = [
        {
            "groups_from_to": (0, num_groups - 1),
            "angular_quadrature": pquad,
            "angle_aggregation_type": "single",
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-6,
            "l_max_its": 30,
            "gmres_restart_interval": 10,
        },
    ]

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=groupsets,
        xs_map=xsecs,
        boundary_conditions=[
            {"name": "Outer", "type": "vacuum"},
            {"name": "Top", "type": "vacuum"},
            {"name": "Bottom", "type": "vacuum"},
        ],
        volumetric_sources=[vol_src],
        options={"save_angular_flux": True}
    )

    # Execute
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Post processing
    bndrys = ["Outer", "Top", "Bottom"]
    leakage = phys.ComputeLeakage(bndrys)
    lkg_outer = leakage['Outer']
    lkg_top = leakage['Top']
    lkg_bottom = leakage['Bottom']

    if rank == 0:
        print(f"Top leakage={(lkg_top).item()}")
        print(f"Bottom leakage={(lkg_bottom).item()}")
        print(f"Side leakage={(lkg_outer).item()}")
        print(f"Total leakage={(lkg_outer + lkg_top + lkg_bottom).item()}")
