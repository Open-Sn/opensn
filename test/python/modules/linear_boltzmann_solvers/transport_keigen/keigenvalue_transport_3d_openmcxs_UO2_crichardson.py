#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 172G KEigenvalue::Solver test using power iteration and OpenMC MGXS library
Test: Final k-eigenvalue: 1.5029618
"""

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
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, NonLinearKEigenSolver

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Cells
    Nx, Ny, Nz = 5, 5, 5
    # Dimensions
    Lx, Ly, Lz = 2.0, 2.0, 2.0

    xmesh = []
    xmin = 0.0
    dx = Lx / Nx
    for i in range(Nx + 1):
        xmesh.append(xmin + i * dx)

    ymesh = []
    ymin = 0.0
    dy = Ly / Ny
    for i in range(Ny + 1):
        ymesh.append(ymin + i * dy)

    zmesh = []
    zmin = 0.0
    dz = Lz / Nz
    for i in range(Nz + 1):
        zmesh.append(zmin + i * dz)

    meshgen = OrthogonalMeshGenerator(node_sets=[xmesh, ymesh, zmesh])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Materials
    xs_uo2 = MultiGroupXS()
    xs_uo2.LoadFromOpenMC("uo2.h5", "set1", 294.0)

    # Solver
    num_groups = 172
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=2,
                    n_azimuthal=4,
                    scattering_order=1
                ),
                "inner_linear_method": "classic_richardson",
                "l_max_its": 500,
                "l_abs_tol": 1.0e-12,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_uo2},
        ],
        scattering_order=1,
        options={
            "boundary_conditions": [
                {"name": "xmin", "type": "reflecting"},
                {"name": "xmax", "type": "reflecting"},
                {"name": "ymin", "type": "reflecting"},
                {"name": "ymax", "type": "reflecting"},
                {"name": "zmin", "type": "reflecting"},
                {"name": "zmax", "type": "reflecting"},
            ],
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": True,
        },
    )

    k_solver = NonLinearKEigenSolver(
        do_problem=phys,
        nl_max_its=500,
        nl_abs_tol=1.0e-8,
    )
    k_solver.Initialize()
    k_solver.Execute()
