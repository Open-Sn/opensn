#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# read a cross section file, apply a scaling factor, and printout some components to the console

# lua function to print entries of a table

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


    function dump(o)
      if type(o) == "table" then
        local s = "{ "
        for k, v in pairs(o) do
          if type(k) ~= "number" then
            k = '"' + k + '"'
          s = s + "[" + k + "] = " + dump(v) + ","
        return s + "} "
      else
        return tostring(o)

    # Create cross sections
    my_xs = []

    my_xs["fuel"] = xs.LoadFromOpenMC(
      "+/+/modules/linear_boltzmann_solvers/transport_keigen/u235_172g.h5",
      "u235",
      294.0
    )

    # print to console
    print("chi\n", dump(my_xs["fuel"].chi))
    print("sigma total\n", dump(my_xs["fuel"].sigma_t))
    chi_before = my_xs["fuel"].chi[1]
    sigt_before = my_xs["fuel"].sigma_t[1]

    # scaling factor
my_xs["fuel"].SetScalingFactor(2.0)

    # print to console again
print("\n After scaling.")
    print("chi\n", dump(my_xs["fuel"].chi))
    print("sigma total\n", dump(my_xs["fuel"].sigma_t))

    chi_after = my_xs["fuel"]["chi"][1]
    sigt_after = my_xs["fuel"]["sigma_t"][1]

log.Log(LOG_0, "chi[1] before. " + chi_before)
log.Log(LOG_0, "chi[1] after . " + chi_after)
log.Log(LOG_0, "sigt[1] before. " + sigt_before)
log.Log(LOG_0, "sigt[1] after . " + sigt_after)
