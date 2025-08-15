#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 2G Infinite Medium Hex import test. Imports EXODUSII.
Uses KEigenvalue::Solver with Power Iteration
Test: Final k-eigenvalue: 0.9293377
"""

# Set and check number of processors

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver

if __name__ == "__main__":

    num_procs = 1
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Setup mesh
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/fuel_hex.e",
    )
    grid = meshgen.Execute()

    # Set Materials (Fuel)
    num_groups = 2
    xs_fuel_g2 = MultiGroupXS()
    xs_fuel_g2.LoadFromOpenSn("xs_fuel_g2.xs")

    # Initialize the LBSProblem
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, 1],
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=8,
                    n_azimuthal=16,
                    scattering_order=1
                ),
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_fuel_g2},
        ],
        scattering_order=1,
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": True,
        },
    )
    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        k_tol=1e-6,
    )
    k_solver.Initialize()
    k_solver.Execute()
