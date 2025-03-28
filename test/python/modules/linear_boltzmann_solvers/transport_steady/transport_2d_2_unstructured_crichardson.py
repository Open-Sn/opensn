#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 2D Transport test with Vacuum and Incident-isotropic BC.
# SDM: PWLD
# Test: Max-value=0.51187 and 1.42458e-03

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationLine, FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":


    num_procs = 4

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Unstructured mesh
    meshgen = MeshGenerator(
      inputs = {
        FromFileMeshGenerator(
          filename = "+/+/+/+/assets/mesh/TriangleMesh2x2Cuts.obj",
        ),
      },
      partitioner = KBAGraphPartitioner(
        nx = 2,
        ny = 2,
        nz = 1,
        xcuts = [ 0.0 ],
        ycuts = [ 0.0 ],
      ),
    )
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume( infx = True, infy = True, infz = True )
    grid.SetUniformBlockID(0)

    num_groups = 168
    xs_3_170 =  MultiGroupXS()
    xs_3_170.LoadFromOpenSn("xs_3_170.xs")

    strength = []
    for g in range(num_groups):
        strength.append(0.)
    mg_src1 = VolumetricSource( block_ids = [ 1 ], group_strength = strength )
    mg_src2 = VolumetricSource( block_ids = [ 2 ], group_strength = strength )

    # Setup Physics
    pquad = GLCProductQuadrature2DXY(8, 32)

    phys = DiscreteOrdinatesSolver(
      mesh = grid,
      num_groups = num_groups,
      groupsets = [
        {
          "groups_from_to": [0, 62],
          "angular_quadrature": pquad,
          "angle_aggregation_num_subsets": 1,
          "inner_linear_method": "classic_richardson",
          "l_abs_tol": 1.0e-6,
          "l_max_its": 1000,
        },
        {
          "groups_from_to": { 63, num_groups - 1 },
          "angular_quadrature": pquad,
          "angle_aggregation_num_subsets": 1,
          "inner_linear_method": "classic_richardson",
          "l_abs_tol": 1.0e-6,
          "l_max_its": 1000,
        },
      ],
      xs_map = [
        { "block_ids": [ 0 ], "xs": xs_3_170 },
      ],
    ]
    bsrc = []
    for g in range(1, num_groups+1):
      bsrc[g] = 0.0
    bsrc[1] = 1.0 / 4.0 / math.pi

    options = {
      "boundary_conditions": [
        {
          "name": "xmin",
          "type": "isotropic",
          "group_strength": bsrc,
        },
      ],
      "scattering_order": 1,
      "volumetric_sources": [ mg_src1, mg_src2 ],
    },
    )


    # Initialize and Execute Solver
    ss_solver = SteadyStateSolver( lbs_solver = phys )

    ss_solver.Initialize()
    ss_solver.Execute()

    # Get field functions
    fflist = GetScalarFieldFunctionList(phys)

    # Volume integrations
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType(OP_MAX)
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[1])

    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()

    if rank == 0:
        print(f"Max-value1={maxval:.5f}")

    # Volume integrations
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType(OP_MAX)
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[10])

    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()

    if rank == 0:
        print(f"Max-value2={maxval:.5e}")
