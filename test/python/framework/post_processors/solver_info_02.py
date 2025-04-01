#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Post-Processor test with lots of post-processors
# Testing table wrapping and getting the value of a post-processor by both
# handle and name

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    # Example Point-Reactor Kinetics solver
    phys0 = PRKSolver(initial_source=0.0)

    # Create many post-processors
    pp = {}
    for k in range(1, 21):
        pp[k] = SolverInfoPostProcessor(
            name="neutron_population" + str(k),
            solver=phys0,
            info={'name': "neutron_population"},
            print_on=[""],
        )

    pp21 = SolverInfoPostProcessor(
        name="neutron_population" + str(20),
        solver=phys0,
        info={'name': "neutron_population"},
        print_on=[""],
    )

    SetPrinterOptions(time_history_limit=5)

    phys0.Initialize()

    for t in range(1, 21):
        phys0.Step()
        time = phys0.GetTimeNew()
        print(t, "{:.3f} {:.5f}".format(time, phys0.GetPopulationNew()), sep="\t")
        phys0.Advance()
        if time > 0.1:
            phys0.SetRho(0.8)

    print("Manual neutron_population1=", f"{pp[1].GetValue():.5f}", sep="\t")
    print("Manual neutron_population1=", f"{pp21.GetValue():.5f}", sep="\t")
