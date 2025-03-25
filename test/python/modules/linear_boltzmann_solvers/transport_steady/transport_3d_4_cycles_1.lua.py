#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 3D Transport test with Vacuum and Incident-isotropic BC.
# SDM: PWLD
# Test: Max-value=3.74343e-04

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
    meshgen = ExtruderMeshGenerator(
      inputs = {
        FromFileMeshGenerator(
          filename = "+/+/+/+/assets/mesh/Square2x2_partition_cyclic3.obj",
        ),
      },
      layers = { { z = 0.4, n = 2 }, { z = 0.8, n = 2 }, { z = 1.2, n = 2 }, { z = 1.6, n = 2 } }, # layers
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
grid.SetBlockIDFromLogicalVolume(vol0, 0, True)

    vol1 =
      RPPLogicalVolume( xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5, infz = True )
grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    num_groups = 21
    xs_graphite =  MultiGroupXS()
    xs_graphite.LoadFromOpenSn("xs_graphite_pure.xs")

    strength = []
    for g in range(1, num_groups+1):
      strength[g] = 0.0

    mg_src1 = VolumetricSource( block_ids = { 1 }, group_strength = strength )
    mg_src2 = VolumetricSource( block_ids = { 2 }, group_strength = strength )

    # Setup Physics
    pquad = GLCProductQuadrature3DXYZ(4, 8)

    lbs_block = {
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, 20 },
          angular_quadrature = pquad,
          #angle_aggregation_type = "single",
          angle_aggregation_num_subsets = 1,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 300,
          gmres_restart_interval = 30,
        },
      },
      xs_map = {
        { block_ids = { 0, 1 }, xs = xs_graphite },
      },
    }
    bsrc = []
    for g in range(1, num_groups+1):
      bsrc[g] = 0.0
    bsrc[1] = 1.0 / 4.0 / math.pi
    lbs_options = {
      boundary_conditions = {
        { name = "zmax", type = "isotropic", group_strength = bsrc },
      },
      scattering_order = 1,
      volumetric_sources = { mg_src1, mg_src2 },
    }
    if reflecting then
      table.insert(lbs_options.boundary_conditions, { name = "zmax", type = "reflecting" )

    phys = DiscreteOrdinatesSolver.Create(lbs_block)
phys.SetOptions(lbs_options)

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
    #    --fieldfunc.SetProperty(slices[k],SLICE_BINORM,{x = 0.0, y = 0.0, y = 1.0)
    #    fieldfunc.Initialize(slices[k])
    #    fieldfunc.Execute(slices[k])
    #    fieldfunc.ExportToPython(slices[k])
    #end

    # Volume integrations
    ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
    curffi = ffi1
curffi.SetOperationType(OP_MAX)
curffi.SetLogicalVolume(vol0)
curffi.AddFieldFunction(fflist[1])

curffi.Initialize()
curffi.Execute()
maxval = curffi.GetValue()

    log.Log(LOG_0, string.format("Max-value1=%.5e", maxval))

    ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
    curffi = ffi1
curffi.SetOperationType(OP_MAX)
curffi.SetLogicalVolume(vol0)
curffi.AddFieldFunction(fflist[20])

curffi.Initialize()
curffi.Execute()
maxval = curffi.GetValue()

    log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))

    # Exports
    if master_export == None then
      fieldfunc.ExportToVTK(fflist[1], "ZPhi3D", "Phi")

      line = fieldfunc.FieldFunctionInterpolationLine.Create()
line.SetInitialPoint({ 0.0, -1.0, 0.5 )
line.SetFinalPoint({ 0.0, 1.0, 0.5 )
line.SetNumberOfPoints(1000)
line.AddFieldFunction(fflist[2])

line.Initialize()
line.Execute()

      fieldfunc.ExportToCSV(line, "Line")

    # Plots
    if location_id == 0 and master_export == None then
#os.system("python ZPFFI00.py")
#--os.system("python ZPFFI11.py")
      #local handle = io.popen("python ZPFFI00.py")
      print("Execution completed")
