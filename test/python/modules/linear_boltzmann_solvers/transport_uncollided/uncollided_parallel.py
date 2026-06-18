#!/usr/bin/env python3

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
run = importlib.import_module("uncollided_parallel_common").run

run(globals(), rank)
