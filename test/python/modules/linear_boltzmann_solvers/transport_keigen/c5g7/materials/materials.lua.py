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


    xss = {}

    for m = 0, 6 do
      xss[tostring(m)] = xs.Create()
    end

    xss["0"] = xs.LoadFromOpenSn("materials/XS_water.xs")
    xss["1"] = xs.LoadFromOpenSn("materials/XS_UO2.xs")
    xss["2"] = xs.LoadFromOpenSn("materials/XS_7pMOX.xs")
    xss["3"] = xs.LoadFromOpenSn("materials/XS_guide_tube.xs")
    xss["4"] = xs.LoadFromOpenSn("materials/XS_4_3pMOX.xs")
    xss["5"] = xs.LoadFromOpenSn("materials/XS_8_7pMOX.xs")
    xss["6"] = xs.LoadFromOpenSn("materials/XS_fission_chamber.xs")

    num_groups = xss["0"].num_groups
    log.Log(LOG_0, "Num groups: " + tostring(num_groups))

    xs_map = {}
    for m = 0, 6 do
      xs_map[m + 1] = { block_ids = { m }, xs = xss[tostring(m)] }
    end
