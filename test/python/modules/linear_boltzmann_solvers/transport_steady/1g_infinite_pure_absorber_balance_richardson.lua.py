#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Infinite, 1-group, pure absorber with balance
# Create Mesh

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

nodes = {}
N = 2
L = 10
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes, nodes } })
grid = meshgen:Execute()

# Set block IDs
grid:SetUniformBlockID(0)

num_groups = 1

# Add cross sections to materials
xs1g = xs.CreateSimpleOneGroup(1.0, 0.0)

strength = {}
strength[1] = 1.0
mg_src = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = strength })

# Angular Quadrature
pquad = aquad.CreateGLCProductQuadrature3DXYZ(4, 8)

# LBS block option
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = pquad,
      inner_linear_method = "petsc_richardson",
      l_abs_tol = 1.0e-9,
      l_max_its = 300,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xs1g },
  },
  options = {
    boundary_conditions = {
      { name = "xmin", type = "reflecting" },
      { name = "xmax", type = "reflecting" },
      { name = "ymin", type = "reflecting" },
      { name = "ymax", type = "reflecting" },
      { name = "zmin", type = "reflecting" },
      { name = "zmax", type = "reflecting" },
    },
    volumetric_sources = { mg_src },
  },
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

# Initialize and execute solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

ss_solver:Initialize()
ss_solver:Execute()

# compute particle balance
phys:ComputeBalance()
