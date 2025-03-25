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
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased
    from pyopensn.fieldfunc import FieldFunctionInterpolationLine, FieldFunctionInterpolationVolume
    from pyopensn.settings import EnableCaliper
    from pyopensn.math import Vector3
    from pyopensn.logvol import RPPLogicalVolume

num_procs = 4
#Unstructured mesh

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
meshgen1 = mesh.MeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "reactor_pin_mesh.obj",
    }),
  },
})
grid = meshgen1:Execute()
# Exports
if master_export == None then
  mesh.ExportToPVTU("ZObjMesh")
end
