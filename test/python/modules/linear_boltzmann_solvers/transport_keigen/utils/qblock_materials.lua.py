#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Create cross sections

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


    xss = []
    xss["0"] = xs.LoadFromOpenSn("xs_water_g2.xs")
    xss["1"] = xs.LoadFromOpenSn("xs_fuel_g2.xs")

    num_groups = xss["0"].num_groups

    # Create materials
    xs_map = []
    for m in range(0, 1+1):
      key = tostring(m)
      xs_map[m + 1] = { block_ids = { m }, xs = xss[key] }
