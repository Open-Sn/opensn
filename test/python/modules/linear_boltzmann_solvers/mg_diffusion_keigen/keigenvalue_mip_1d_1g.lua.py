#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 1D 1G KEigenvalue::Solver test using power iteration
# Test: Final k-eigenvalue: 0.9995433

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

# NOTE: For command line inputs, specify as:
#       variable=[[argument]]

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

MPIBarrier()

# ##################################################
# ##### Parameters #####
# ##################################################

# Mesh variables
if L == None then
  L = 100.0
end
if n_cells == None then
  n_cells = 50
end

# Transport angle information
if n_angles == None then
  n_angles = 32
end
if scat_order == None then
  scat_order = 0
end

# k-eigenvalue iteration parameters
if kes_max_iterations == None then
  kes_max_iterations = 5000
end
if kes_tolerance == None then
  kes_tolerance = 1e-8
end

# Source iteration parameters
if si_max_iterations == None then
  si_max_iterations = 500
end
if si_tolerance == None then
  si_tolerance = 1e-8
end

# Delayed neutrons
if use_precursors == None then
  use_precursors = True
end

# ##################################################
# ##### Run problem #####
# ##################################################

# Setup mesh
nodes = {}
dx = L / n_cells
for i = 0, n_cells do
  nodes[i + 1] = i * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
grid = meshgen1:Execute()

# Set block IDs
grid:SetUniformBlockID(0)

xs_simple_fissile = xs.LoadFromOpenSn("+/transport_keigen/simple_fissile.xs")

# Setup Physics
num_groups = 1
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = aquad.CreateGLProductQuadrature1DSlab(n_angles),
      inner_linear_method = "petsc_gmres",
      l_max_its = si_max_iterations,
      l_abs_tol = si_tolerance,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xs_simple_fissile },
  },
}

lbs_options = {
  scattering_order = scat_order,

  use_precursors = use_precursors,

  verbose_inner_iterations = False,
  verbose_outer_iterations = True,
}

phys = lbs.DiffusionDFEMSolver.Create(lbs_block)
phys:SetOptions(lbs_options)

k_solver0 = lbs.NonLinearKEigen.Create({
  lbs_solver = phys,
  nl_max_its = kes_max_iterations,
  nl_abs_tol = kes_tolerance,
})
k_solver0:Initialize()
k_solver0:Execute()

# Get field functions
# Line plot
# Volume integrations
# Exports
# Plots
