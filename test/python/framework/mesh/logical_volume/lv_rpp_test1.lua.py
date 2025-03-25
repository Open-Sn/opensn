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


    lv1 = logvol.RPPLogicalVolume( xmin = -12.0, xmax = 12.0 )

print("lv1 test 1.", lv1.Inside({ x = -13.0, y = 0.1, z = 0.1 ))
print("lv1 test 2.", lv1.Inside({ x = -11.0, y = 0.1, z = 0.1 ))
print("lv1 test 3.", lv1.Inside({ x = -11.0, y = -1.0, z = -1.0 ))

    lv2 = logvol.RPPLogicalVolume( xmin = -12.0, xmax = 12.0, infy = True, infz = True )
    print()
print("lv2 test 1.", lv2.Inside({ x = -13.0, y = 0.1, z = 0.1 ))
print("lv2 test 2.", lv2.Inside({ x = -11.0, y = 0.1, z = 0.1 ))
print("lv2 test 3.", lv2.Inside({ x = -11.0, y = -1.0, z = -1.0 ))
