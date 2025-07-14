#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D Transport test with Vacuum and Incident-isotropic BC.
SDM: PWLD
Test: Max-value=3.74343e-04
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume, FieldFunctionInterpolationLine
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3

if __name__ == "__main__":

    num_procs = 4

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    meshgen = ExtruderMeshGenerator(
        inputs=[
            FromFileMeshGenerator(
                filename="../../../../assets/mesh/Square2x2_partition_cyclic3.obj"
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
    vol1 = RPPLogicalVolume(xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    # Cross sections
    num_groups = 21
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("xs_graphite_pure.xs")

    # Source
    strength = [0.0 for _ in range(num_groups)]
    mg_src1 = VolumetricSource(block_ids=[1], group_strength=strength)
    mg_src2 = VolumetricSource(block_ids=[2], group_strength=strength)

    # Setup Physics
    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=1)

    bsrc = [0.0 for _ in range(num_groups)]
    bsrc[0] = 1.0

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, 20],
                "angular_quadrature": pquad,
                # "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0, 1], "xs": xs_graphite},
        ],
        scattering_order=1,
        options={
            "boundary_conditions": [
                {"name": "zmax", "type": "isotropic", "group_strength": bsrc},
            ],
            "volumetric_sources": [mg_src1, mg_src2],
        },
    )

    # Initialize and Execute Solver
    ss_solver = SteadyStateSolver(lbs_problem=phys)
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
    curffi.AddFieldFunction(fflist[19])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")

    line = FieldFunctionInterpolationLine()
    line.SetInitialPoint(Vector3(0.0, -1.0, 0.5))
    line.SetFinalPoint(Vector3(0.0, 1.0, 0.5))
    line.SetNumberOfPoints(1000)
    line.AddFieldFunction(fflist[1])
    line.Initialize()
    line.Execute()
