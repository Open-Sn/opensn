#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Setup mesh

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

meshgen1 = mesh.ExtruderMeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "+/+/+/+/+/assets/mesh/TriangleMesh2x2.obj",
    }),
  },
  layers = { { z = 0.4, n = 2 }, { z = 0.8, n = 2 }, { z = 1.2, n = 2 }, { z = 1.6, n = 2 } }, # layers
})
grid = meshgen1:Execute()
