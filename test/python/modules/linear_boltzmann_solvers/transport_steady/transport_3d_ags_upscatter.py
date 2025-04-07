#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D Transport test with multiple groupsets and groupset-to-groupset upscattering
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
    from pyopensn.mesh import FromFileMeshGenerator, ExtruderMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased
    from pyopensn.fieldfunc import FieldFunctionInterpolationLine, FieldFunctionInterpolationVolume
    from pyopensn.settings import EnableCaliper
    from pyopensn.math import Vector3
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    num_procs = 4

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    meshgen = ExtruderMeshGenerator(
        inputs=[
            FromFileMeshGenerator(
                filename="../../../../assets/mesh/TriangleMesh2x2Cuts.obj"
            )
        ],
        layers=[{"z": 0.4, "n": 2},
                {"z": 0.8, "n": 2},
                {"z": 1.2, "n": 2},
                {"z": 1.6, "n": 2},
                ],  # layers
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            xcuts=[0.0],
            ycuts=[0.0], ),
    )
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)

    # Cross sections
    num_groups = 3
    xs_upscatter = MultiGroupXS()
    xs_upscatter.LoadFromOpenSn("simple_upscatter.xs")

    # Source
    strength = [1.0 for _ in range(num_groups)]
    mg_src = VolumetricSource(block_ids=[0], group_strength=strength)

    # Setup Physics
    pquad = GLCProductQuadrature3DXYZ(4, 8)

    phys = DiscreteOrdinatesSolver(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
            {
                "groups_from_to": [1, 1],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
            {
                "groups_from_to": [2, 2],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_upscatter},
        ],
        options={
            "scattering_order": 0,
            "verbose_ags_iterations": True,
            "max_ags_iterations": 30,
            "ags_tolerance": 1.0e-6,
            "volumetric_sources": [mg_src],
        },
    )
    ss_solver = SteadyStateSolver(lbs_solver=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Get field functions
    fflist = phys.GetScalarFieldFunctionList()

    # Volume integrations
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value1={maxval:.5e}")

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[1])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[2])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value3={maxval:.5e}")
