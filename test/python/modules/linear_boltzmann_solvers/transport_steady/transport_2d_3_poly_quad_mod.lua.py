#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 2D Transport test with Vacuum and Incident-isotropic BC. Quadrature optimized for polar symmetry.
# SDM: PWLD
# Test: Max-value=0.50758 and 2.52527e-04

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
    from pyopensn.aquad import GLProductQuadrature1DSlab
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
    meshgen = MeshGenerator(
      inputs = {
        FromFileMeshGenerator(
          filename = "+/+/+/+/assets/mesh/SquareMesh2x2QuadsBlock.obj",
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
    grid:SetUniformBlockID(0)

    num_groups = 168
    xs_3_170 =  MultiGroupXS()
    xs_3_170.LoadFromOpenSn("xs_3_170.xs")

    strength = []
    for g in range(1, num_groups+1):
      strength[g] = 0.0
    #src[1] = 1.0
    mg_src1 = VolumetricSource( block_ids = { 1 }, group_strength = strength )
    mg_src2 = VolumetricSource( block_ids = { 2 }, group_strength = strength )

    # Setup Physics
    pquad = GLCProductQuadrature2DXY(2, 8)

    lbs_block = {
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, 62 },
          angular_quadrature = pquad,
          angle_aggregation_num_subsets = 1,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 300,
          gmres_restart_interval = 100,
        },
        {
          groups_from_to = { 63, num_groups - 1 },
          angular_quadrature = pquad,
          angle_aggregation_num_subsets = 1,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 300,
          gmres_restart_interval = 100,
        },
      },
      xs_map = [
        { "block_ids": [ 0 ], "xs": xs_3_170 },
      ],
    }
    bsrc = []
    for g in range(1, num_groups+1):
      bsrc[g] = 0.0
    bsrc[1] = 1.0 / 4.0 / math.pi

    lbs_options = {
      boundary_conditions = {
        {
          name = "xmin",
          type = "isotropic",
          group_strength = bsrc,
        },
      },
      scattering_order = 1,
      max_ags_iterations = 1,
      volumetric_sources = { mg_src1, mg_src2 },
    }

    phys = DiscreteOrdinatesSolver.Create(lbs_block)
phys.SetOptions(lbs_options)

    # Initialize and Execute Solver
    ss_solver = SteadyStateSolver( lbs_solver = phys )

ss_solver.Initialize()
ss_solver.Execute()

    # Get field functions
    fflist = GetScalarFieldFunctionList(phys)

    # Volume integrations
    vol0 = RPPLogicalVolume( infx = True, infy = True, infz = True )

    ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
    curffi = ffi1
curffi.SetOperationType(OP_MAX)
curffi.SetLogicalVolume(vol0)
curffi.AddFieldFunction(fflist[1])

curffi.Initialize()
curffi.Execute()
maxval = curffi.GetValue()

    log.Log(LOG_0, string.format("Max-value1=%.5f", maxval))

    # Volume integrations
    ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
    curffi = ffi1
curffi.SetOperationType(OP_MAX)
curffi.SetLogicalVolume(vol0)
curffi.AddFieldFunction(fflist[160])

curffi.Initialize()
curffi.Execute()
maxval = curffi.GetValue()

    log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))
