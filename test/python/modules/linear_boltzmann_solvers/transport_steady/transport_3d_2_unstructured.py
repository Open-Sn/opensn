#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D PWLD unstructured mesh ransport test with vacuum and incident-isotropic boundary conditions
Test: Max-value=5.41465e-01 and 3.78243e-04
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
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Setup mesh
    meshgen = ExtruderMeshGenerator(
        inputs=[
            FromFileMeshGenerator(
                filename="../../../../assets/mesh/TriangleMesh2x2Cuts.obj"
            )
        ],
        layers=[
            {"z": 0.4, "n": 2},
            {"z": 0.8, "n": 2},
            {"z": 1.2, "n": 2},
            {"z": 1.6, "n": 2},
        ],
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            xcuts=[0.0],
            ycuts=[0.0],
        )
    )
    grid = meshgen.Execute()

    # Set block IDs using logical volumes
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)
    vol1 = RPPLogicalVolume(xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    # Energy groups and cross-section data.
    num_groups = 21
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("xs_graphite_pure.xs")

    # Define sources
    strength = [0.0 for _ in range(num_groups)]
    mg_src1 = VolumetricSource(block_ids=[1], group_strength=strength)
    mg_src2 = VolumetricSource(block_ids=[2], group_strength=strength)

    # Angular quadrature
    pquad = GLCProductQuadrature3DXYZ(4, 8)

    # Set up the boundary source.
    bsrc = [0.0 for _ in range(num_groups)]
    bsrc[0] = 1.0 / 4.0 / math.pi

    # Create and configure the discrete ordinates solver
    phys = DiscreteOrdinatesSolver(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, 20),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0, 1], "xs": xs_graphite},
        ],
        options={
            "boundary_conditions": [
                {"name": "zmin", "type": "isotropic", "group_strength": bsrc},
            ],
            "scattering_order": 1,
            "save_angular_flux": True,
            "volumetric_sources": [mg_src1, mg_src2],
        }
    )
    ss_solver = SteadyStateSolver(lbs_solver=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Retrieve the scalar field function list.
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)

    # Volume integration for the first field function
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0][0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value1={maxval:.5e}")

    # Volume integration for the twentieth field function
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[19][0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")
