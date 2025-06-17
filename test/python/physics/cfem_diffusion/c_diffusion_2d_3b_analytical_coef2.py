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
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased
    from pyopensn.fieldfunc import FieldFunctionInterpolationLine, FieldFunctionInterpolationVolume
    from pyopensn.settings import EnableCaliper
    from pyopensn.math import Vector3
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    # Setup mesh
    nodes = []
    N = 40
    L = 1
    xmin = 0
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    # governing law: -(u_xx + u_yy) = q, on domain [0,1]x[0,1]
    # when the exact solution is chosen u(x,y) = sin(pi.x) * sin(pi.y)
    # this automatically gives:
    #    boundary = zero-Dirichlet on all 4 sides
    #    volumetric source term: q(,x) = 2*pi*pi * sin(pi.x) * sin(pi.y)
    # the factor 2 is the dim of the problem

    def D_coef(i, pt):
        return 1.0

    def Q_ext(i, pt):
        return 2. * math.pi * math.pi * math.sin(math.pi * pt.x) * math.sin(math.pi * pt.y)

    def Sigma_a(i, pt):
        return 0.0

    # Setboundary IDs
    # xmin,xmax,ymin,ymax,zmin,zmax
    e_vol = RPPLogicalVolume(xmin=0.99999, xmax=1000.0, infy=True, infz=True)
    w_vol = RPPLogicalVolume(xmin=-1000.0, xmax=-0.99999, infy=True, infz=True)
    n_vol = RPPLogicalVolume(ymin=0.99999, ymax=1000.0, infx=True, infz=True)
    s_vol = RPPLogicalVolume(ymin=-1000.0, ymax=-0.99999, infx=True, infz=True)

    e_bndry = "0"
    w_bndry = "1"
    n_bndry = "2"
    s_bndry = "3"

    grid.SetBoundaryIDFromLogicalVolume(e_vol, e_bndry, True)
    grid.SetBoundaryIDFromLogicalVolume(w_vol, w_bndry, True)
    grid.SetBoundaryIDFromLogicalVolume(n_vol, n_bndry, True)
    grid.SetBoundaryIDFromLogicalVolume(s_vol, s_bndry, True)

    d_coef_fn = ScalarSpatialMaterialFunction(D_coef)
    Q_ext_fn = ScalarSpatialMaterialFunction(Q_ext)
    Sigma_a_fn = ScalarSpatialMaterialFunction(Sigma_a)

    # CFEM solver
    phys = CFEMDiffusionSolver(
        name="CFEMDiffusionSolver",
        mesh=grid,
        residual_tolerance=1e-8,
    )
    phys.SetOptions(boundary_conditions=[
        {"boundary": e_bndry, "type": "dirichlet", "coeffs": [0.0]},
        {"boundary": n_bndry, "type": "dirichlet", "coeffs": [0.0]},
        {"boundary": s_bndry, "type": "dirichlet", "coeffs": [0.0]},
        {"boundary": w_bndry, "type": "dirichlet", "coeffs": [0.0]}
    ])
    phys.SetDCoefFunction(d_coef_fn)
    phys.SetQExtFunction(Q_ext_fn)
    phys.SetSigmaAFunction(Sigma_a_fn)
    phys.Initialize()
    phys.Execute()

    # Get field functions
    fflist = phys.GetFieldFunctions()

    # PostProcessors
    maxval = AggregateNodalValuePostProcessor(
        name="maxval",
        field_function=fflist[0],
        operation="max",
    )
    maxval.Execute()
