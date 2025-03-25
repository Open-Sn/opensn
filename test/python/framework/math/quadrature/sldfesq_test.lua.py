#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

local function split_string(input, separator)
  local result = {}
  for value in input:gmatch("[^" + separator + "]+") do
    table.insert(result, value)
  end
  return result
end

local function getWeightSum(quad_points)
  log.Log(LOG_0, string.format("\nReading : %s", quad_points))
  local weight_sum = 0.0
  local file = io.open(quad_points, "r")
  if file then
    for line in file:lines() do
      local values = split_string(line, " ")
      local float_values = {}
      weight_sum = weight_sum + tonumber(values[4])
    end
    file:close()
  else
    log.Log(LOG_0ERROR, string.format("Error: Could not open file %s", quad_points))
  end
  return weight_sum
end

# Qudrature-1 : Test creation of SLDFESQ with initial refinement level of 0
cquad1 = aquad.CreateSLDFESQuadrature(0)
aquad.PrintQuadratureToFile(cquad1, "TestQuad1_")
local quad1_sum = getWeightSum("TestQuad1_points.txt")
log.Log(LOG_0, string.format("Weight-Sum-1=%.3e\n\n", quad1_sum / (4 * math.pi)))

# Qudrature-2 : Test local refinement of SLDFESQ with initial refinement level of 1
cquad2 = aquad.CreateSLDFESQuadrature(2)
aquad.LocallyRefineSLDFESQ(cquad2, { x = 0.25, y = -0.85, z = 1.0 }, 30.0 * math.pi / 180, False)
aquad.PrintQuadratureToFile(cquad2, "TestQuad2_")
local quad2_sum = getWeightSum("TestQuad2_points.txt")
log.Log(LOG_0, string.format("Weight-Sum-2=%.3e", quad2_sum / (4 * math.pi)))

os.execute("rm TestQuad1* TestQuad2*")
