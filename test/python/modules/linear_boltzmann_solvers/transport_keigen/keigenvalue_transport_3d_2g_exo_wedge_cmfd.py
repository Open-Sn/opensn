#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 2G EXODUS wedge k-eigenvalue test using identity CMFD acceleration.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration

if __name__ == "__main__":

    num_procs = 1
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/fuel_wedge.e",
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()

    num_groups = 2
    xs_fuel_g2 = MultiGroupXS()
    xs_fuel_g2.LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")

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
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )
    cmfd = CMFDAcceleration(
        problem=phys,
        coarse_mesh="identity",
        relaxation=1.0,
        l_abs_tol=1.0e-10,
        max_iters=100,
        verbose=False,
        petsc_options="-CMFDAccelerationksp_type gmres -CMFDAccelerationpc_type jacobi",
    )
    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        acceleration=cmfd,
        k_tol=1e-6,
    )
    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()
    if rank == 0:
        print(f"Python k-eigenvalue: {k}")

    if abs(k - 0.9293377) > 5.0e-6:
        raise RuntimeError(f"Expected k near 0.9293377, got {k}")
