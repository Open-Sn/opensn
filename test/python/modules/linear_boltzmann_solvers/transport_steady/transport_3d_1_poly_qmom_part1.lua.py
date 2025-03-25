#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 3D Transport test with Vacuum BCs and a material source writing source moments.
# SDM: PWLD
# Test: Max-value=1.08320e-01 and 0.0

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

num_procs = 4

if size != num_procs:
    sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

# Setup mesh
meshgen1 = mesh.ExtruderMeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "+/+/+/+/assets/mesh/SquareMesh2x2Quads.obj",
    }),
  },
  layers = { { z = 0.4, n = 2 }, { z = 0.8, n = 2 }, { z = 1.2, n = 2 }, { z = 1.6, n = 2 } }, # layers
  partitioner = mesh.KBAGraphPartitioner.Create({
    nx = 2,
    ny = 2,
    xcuts = { 0.0 },
    ycuts = { 0.0 },
  }),
})
grid = meshgen1:Execute()

# Set block IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = True, infy = True, infz = True })
grid:SetBlockIDFromLogicalVolume(vol0, 0, True)

vol1 = logvol.RPPLogicalVolume.Create({
  xmin = -0.5 / 8,
  xmax = 0.5 / 8,
  ymin = -0.5 / 8,
  ymax = 0.5 / 8,
  infz = True,
})
grid:SetBlockIDFromLogicalVolume(vol1, 1, True)

num_groups = 21
xs_graphite = xs.LoadFromOpenSn("xs_graphite_pure.xs")

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
mg_src0 = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = strength })
strength[1] = 1.0
mg_src1 = lbs.VolumetricSource.Create({ block_ids = { 1 }, group_strength = strength })

# Setup Physics
pquad0 = aquad.CreateGLCProductQuadrature3DXYZ(4, 8)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 20 },
      angular_quadrature = pquad0,
      #angle_aggregation_type = "single",
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_richardson",
      l_abs_tol = 1.0e-6,
      l_max_its = 2,
      gmres_restart_interval = 100,
    },
  },
  xs_map = {
    { block_ids = { 0, 1 }, xs = xs_graphite },
  },
}

lbs_options = {
  scattering_order = 1,
  volumetric_sources = { mg_src0, mg_src1 },
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

# Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

lbs.CreateAndWriteSourceMoments(phys1, "Qmoms")

# Get field functions
fflist = lbs.GetScalarFieldFunctionList(phys1)

################################################ Exports

if master_export == None then
  fieldfunc.ExportToVTKMulti(fflist, "ZPhi3D")
end

# Plots
if location_id == 0 and master_export == None then
  #os.execute("python ZPFFI00.py")
  #--os.execute("python ZPFFI11.py")
  #local handle = io.popen("python ZPFFI00.py")
  print("Execution completed")
end
