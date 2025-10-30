#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D Transport test with distributedt-mesh + 2D ortho mesh + extruded mesh.
SDM: PWLD
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator, DistributedMeshGenerator
    from pyopensn.mesh import ExtruderMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Cells
    div = 8
    Nx = math.floor(128 / div)
    Ny = math.floor(128 / div)
    Nz = math.floor(256 / div)

    # Dimensions
    Lx = 10.0
    Ly = 10.0
    Lz = 10.0

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

    meshgen = DistributedMeshGenerator(
        inputs=[
            OrthogonalMeshGenerator(node_sets=[xmesh, ymesh]),
            ExtruderMeshGenerator(
                layers=[{"z": Lz, "n": Nz}],  # layers
            ),
        ],
    )

    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)

    # Cross sections
    num_groups = 21
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("xs_graphite_pure.xs")

    # Source
    strength = [0.0 for _ in range(num_groups)]
    mg_src = VolumetricSource(block_ids=[1], group_strength=strength)

    # Setup Physics
    pquad = GLCProductQuadrature3DXYZ(n_polar=8, n_azimuthal=8, scattering_order=1)

    bsrc = [0.0 for _ in range(num_groups)]
    bsrc[0] = 1.0

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, 20],
                "angular_quadrature": pquad,
                "angle_aggregation_type": "polar",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_graphite},
        ],
        scattering_order=1,
        volumetric_sources=[mg_src],
        boundary_conditions=[
            {"name": "xmin", "type": "isotropic", "group_strength": bsrc},
        ],
        options={
            "save_angular_flux": True,
        },
    )

    # Initialize and Execute Solver
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Get field functions
    fflist = phys.GetScalarFieldFunctionList()

    ffi1 = FieldFunctionInterpolationVolume()
    ffi1.SetOperationType("max")
    ffi1.SetLogicalVolume(vol0)
    ffi1.AddFieldFunction(fflist[0])
    ffi1.Initialize()
    ffi1.Execute()
    maxval = ffi1.GetValue()
    if rank == 0:
        print(f"Max-value-0={maxval:.5e}")

    ffi1 = FieldFunctionInterpolationVolume()
    ffi1.SetOperationType("max")
    ffi1.SetLogicalVolume(vol0)
    ffi1.AddFieldFunction(fflist[19])
    ffi1.Initialize()
    ffi1.Execute()
    maxval = ffi1.GetValue()
    if rank == 0:
        print(f"Max-value-19={maxval:.5e}")
