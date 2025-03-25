#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 3D Transport test with Vacuum and Incident-isotropic BC.
# SDM: PWLD
# Test: Max-value=0.49903 and 7.18243e-4

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

# Check num_procs
if check_num_procs == None and number_of_processes ~= num_procs then
  log.Log(
    LOG_0ERROR,
    "Incorrect amount of processors. "
      + "Expected "
      + tostring(num_procs)
      + ". Pass check_num_procs=False to override if possible."
  )
  os.exit(False)
end

# Setup mesh
Nxy = 32
nodesxy = {}
dxy = 2 / Nxy
dz = 1.6 / 8
for i = 0, Nxy do
  nodesxy[i + 1] = -1.0 + i * dxy
end
nodesz = {}
for k = 0, 8 do
  nodesz[k + 1] = 0.0 + k * dz
end

meshgen = mesh.MeshGenerator.Create({
  inputs = {
    mesh.OrthogonalMeshGenerator.Create({
      node_sets = { nodesxy, nodesxy, nodesz },
    }),
  },
  partitioner = mesh.PETScGraphPartitioner.Create({ type = "ptscotch" }),
})
grid = meshgen:Execute()

# Set block IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = True, infy = True, infz = True })
grid:SetBlockIDFromLogicalVolume(vol0, 0, True)

vol1 =
  logvol.RPPLogicalVolume.Create({ xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5, infz = True })
grid:SetBlockIDFromLogicalVolume(vol1, 1, True)

num_groups = 21
xs_graphite = xs.LoadFromOpenSn("xs_graphite_pure.xs")

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
mg_src1 = lbs.VolumetricSource.Create({ block_ids = { 1 }, group_strength = strength })
mg_src2 = lbs.VolumetricSource.Create({ block_ids = { 2 }, group_strength = strength })

# Setup Physics
pquad0 = aquad.CreateGLCProductQuadrature3DXYZ(4, 8)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 20 },
      angular_quadrature = pquad0,
      angle_aggregation_type = "single",
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
  xs_map = {
    { block_ids = { 0, 1 }, xs = xs_graphite },
  },
}
bsrc = {}
for g = 1, num_groups do
  bsrc[g] = 0.0
end
bsrc[1] = 1.0 / 4.0 / math.pi
lbs_options = {
  boundary_conditions = {
    { name = "zmin", type = "isotropic", group_strength = bsrc },
  },
  scattering_order = 1,
  volumetric_sources = { mg_src1, mg_src2 },
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

# Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

# Get field functions
fflist = lbs.GetScalarFieldFunctionList(phys1)

# Volume integrations
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[1])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value1=%.5e", maxval))

ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[20])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))
