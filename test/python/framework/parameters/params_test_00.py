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

if __name__ == "__main__":


    block = {
      enabled = True,
      it_method = "gmres",
      nl_abs_tol = 1.0e-12,
      nl_max_its = 33,
      sub1 = {
        ax_method = 2,
        "l_abs_tol": 1.0e-2,
      },
      sub2 = {
        ax_method = 3,
        "l_abs_tol": 1.0e-3,
        blocks = { 99, 98, 97 },
        cblocks = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } },
      },
    }

    unit_tests.ParameterBlock_Test00(block)
