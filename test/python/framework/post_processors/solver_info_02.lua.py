#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Post-Processor test with lots of post-processors
# Testing table wrapping and getting the value of a post-processor by both
# handle and name

# Example Point-Reactor Kinetics solver

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


    phys0 = prk.PRKSolver( initial_source = 0.0 )

    pp = []
    for k in range(1, 20+1):
      pp[k] = post.SolverInfoPostProcessor(
        name = "neutron_population" + tostring(k),
        solver = phys0,
        info = { name = "neutron_population" },
        print_on = { "" },
      )
    pp21 = post.SolverInfoPostProcessor(
      name = "neutron_population" + tostring(21),
      solver = phys0,
      info = { name = "neutron_population" },
      print_on = { "" },
    )

    post.SetPrinterOptions({
      time_history_limit = 5,
    )

    solver.Initialize(phys0)

    for t in range(1, 20+1):
      solver.Step(phys0)
time = phys0.TimeNew()
print(t, string.format("%.3f %.5f", time, phys0.PopulationNew()))

      solver.Advance(phys0)
      if time > 0.1 then
phys0.SetRho(0.8)

print("Manual neutron_population1=", string.format("%.5f", pp[1].GetValue()))
print("Manual neutron_population1=", string.format("%.5f", pp21.GetValue()))
