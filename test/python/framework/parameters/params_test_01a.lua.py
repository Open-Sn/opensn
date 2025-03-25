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


    sub_obj = {
      num_groups = 2,
    }

    #Optional parameter "limiter_type". Should create a deprecation warning
    unit_testsB.TestObject.Create({
      solver_type = "B",
      coupled_field = "T",
      sub_obj1 = sub_obj,
      limiter_type = 2,
    })

    #Optional parameter "scheme". Should create a deprecation error.
    unit_testsB.TestObject.Create({
      solver_type = "B",
      coupled_field = "T",
      sub_obj1 = sub_obj,
      scheme = "Snotty",
    })
