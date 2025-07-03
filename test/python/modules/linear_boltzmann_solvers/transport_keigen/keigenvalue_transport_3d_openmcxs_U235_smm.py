#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 84G keigenvalue test using OpenMC MGXS cross-sections, power iteration, and SMM
Test: Final k-eigenvalue: 2.280431
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
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import SMMAcceleration

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Mesh setup
    Nx = 5
    Lx = 2.0
    xmin = 0.0
    dx = Lx / Nx
    xmesh = [xmin + k * dx for k in range(Nx + 1)]

    Ny = 5
    Ly = 2.0
    ymin = 0.0
    dy = Ly / Ny
    ymesh = [ymin + k * dy for k in range(Ny + 1)]

    Nz = 5
    Lz = 2.0
    zmin = 0.0
    dz = Lz / Nz
    zmesh = [zmin + k * dz for k in range(Nz + 1)]

    meshgen = OrthogonalMeshGenerator(node_sets=[xmesh, ymesh, zmesh])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Load cross section data using the OpenMC MGXS library
    num_groups = 84
    xs_u235 = MultiGroupXS()
    xs_u235.LoadFromOpenMC("u235_84g.h5", "set1", 294.0)

    # Solver Setup
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
                "l_max_its": 2,
                "l_abs_tol": 1.0e-12,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_u235},
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
        }
    )
    smm = SMMAcceleration(
        lbs_problem=phys,
        verbose=True,
        pi_k_tol=1.0e-8,
        pi_max_its=30,
        sdm="pwld"
    )
    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        lbs_acceleration=smm,
        k_tol=1.0e-8,
    )
    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()
    print(f"Python k-eigenvalue: {k}")
