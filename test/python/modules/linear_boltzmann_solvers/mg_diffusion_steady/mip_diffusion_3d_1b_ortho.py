#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 3D Transport test with Vacuum and Incident-isotropic BC.
# SDM: PWLD
# Test: Max-value=5.28310e-01 and 8.04576e-04

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
    if reflecting == None then
      reflecting = True

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    nodes = []
    N = 10
    L = 5
    xmin = -L / 2
    dx = L / N
    for i in range(N+1):
      nodes.append(xmin + i * dx)
    znodes = []
    for i in range(1, (N / 2 + 1)+1):
      znodes.append(xmin + i * dx)
    if reflecting then
      meshgen = OrthogonalMeshGenerator( node_sets = { nodes, nodes, znodes } )
    else
      meshgen = OrthogonalMeshGenerator(node_sets = [nodes, nodes, nodes])
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume( infx = True, infy = True, infz = True )
    grid:SetUniformBlockID(0)

    num_groups = 21
    xs_graphite =  MultiGroupXS()
    xs_graphite.LoadFromOpenSn("+/transport_steady/xs_graphite_pure.xs")

    strength = []
    for g in range(1, num_groups+1):
      strength[g] = 0.0
    strength[1] = 1.0
    mg_src = VolumetricSource( block_ids = [ 0 ], group_strength = strength )

    # Setup Physics
    pquad = GLCProductQuadrature2DXY(4, 8)

    phys = DiscreteOrdinatesSolver(
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = [0, 20],
          angular_quadrature = pquad,
          angle_aggregation_type = "single",
          angle_aggregation_num_subsets = 1,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 300,
          gmres_restart_interval = 100,
        },
      },
      xs_map = [
        { "block_ids": [ 0 ], "xs": xs_graphite },
      ],
    ]

    lbs_options = [
      "scattering_order": 1,
      "volumetric_sources": [ mg_src ],
    ]
    if reflecting then
      lbs_options."boundary_conditions": [ { "name": "zmax", "type": "reflecting" } ]

    phys = DiffusionDFEMSolver.Create(lbs_block)

    # Initialize and Execute Solver
    ss_solver = SteadyStateSolver( lbs_solver = phys )

    ss_solver.Initialize()
    ss_solver.Execute()

    # Get field functions
    fflist = GetScalarFieldFunctionList(phys)

    # Slice plot
    #slices = []
    #for k in range(1, count+1):
    #    slices[k] = fieldfunc.FFInterpolationCreate(SLICE)
    #    fieldfunc.SetProperty(slices[k],SLICE_POINT,{x = 0.0, y = 0.0, z = 0.8001)
    #    fieldfunc.SetProperty(slices[k],ADD_FIELDFUNCTION,fflist[k])
    #    --fieldfunc.SetProperty(slices[k],SLICE_TANGENT,{x = 0.393, y = 1.0-0.393, z = 0)
    #    --fieldfunc.SetProperty(slices[k],SLICE_NORMAL,{x = -(1.0-0.393), y = -0.393, z = 0.0)
    #    --fieldfunc.SetProperty(slices[k],SLICE_BINORM,{x = 0.0, y = 0.0, z = 1.0)
    #    fieldfunc.Initialize(slices[k])
    #    fieldfunc.Execute(slices[k])
    #    fieldfunc.ExportToPython(slices[k])
    #end

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
    print(f"Max-value1={maxval:.5e}")

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType(OP_MAX)
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[20])

    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()

    if rank == 0:
    print(f"Max-value2={maxval:.5e}")

    # Exports
    if master_export == None then
      if reflecting then
        fieldfunc.ExportToVTKMulti(fflist, "ZPhi3DReflected")
      else
        fieldfunc.ExportToVTKMulti(fflist, "ZPhi3D")

    # Plots
    if location_id == 0 and master_export == None then
#os.system("python ZPFFI00.py")
#--os.system("python ZPFFI11.py")
      #local handle = io.popen("python ZPFFI00.py")
      print("Execution completed")

    # DO
    #[0]  Max-value1=2.52092e+00
    #[0]  Max-value2=5.79100e-03

    # MGDiffusion
    #[0]  Max-value1=2.74873e-01
    #[0]  Max-value2=9.47508e-05
