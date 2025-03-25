#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 2D test to compute and display volume per material ID.

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

num_procs = 3

if size != num_procs:
    sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

# ############################################### Setup mesh
nodes = {}
N = 20
L = 5.
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen1:Execute()

# ############################################### Set Material IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = True, infy = True, infz = True })
grid:SetBlockIDFromLogicalVolume(vol0, 0, True)
vol1 = logvol.RPPLogicalVolume.Create({ xmin = -1000.0, xmax = L / N, infy = True, infz = True })
grid:SetBlockIDFromLogicalVolume(vol1, 1, True)

grid:ComputeVolumePerBlockID()
