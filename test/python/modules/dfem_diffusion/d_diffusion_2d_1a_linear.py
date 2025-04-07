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

    # Setup mesh
    nodes = []
    N = 10
    L = 2
    xmin = -L / 2
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    D = [1.0]
    Q = [0.0]
    XSa = [0.0]

    def D_coef(i, pt):
        return D[i]

    def Q_ext(i, pt):
        return Q[i]

    def Sigma_a(i, pt):
        return XSa[i]

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

    # DFEM solver
    phys = DFEMDiffusionSolver(
        name="CFEMDiffusionSolver",
        mesh=grid,
        residual_tolerance=1e-8,
    )
    phys.SetOptions(boundary_conditions=[
        {"boundary": e_bndry, "type": "robin", "coeffs": [0.25, 0.5, 0.0]},
        {"boundary": n_bndry, "type": "reflecting"},
        {"boundary": s_bndry, "type": "reflecting"},
        {"boundary": w_bndry, "type": "robin", "coeffs": [0.25, 0.5, 1.0]}
    ])
    phys.SetDCoefFunction(d_coef_fn)
    phys.SetQExtFunction(Q_ext_fn)
    phys.SetSigmaAFunction(Sigma_a_fn)
    phys.Initialize()
    phys.Execute()

    # Get field functions
    fflist = phys.GetFieldFunctions()

    # Line plot
    cline = FieldFunctionInterpolationLine()
    cline.SetInitialPoint(Vector3(-L / 2, 0.0, 0.0))
    cline.SetFinalPoint(Vector3(L / 2, 0.0, 0.0))
    cline.SetNumberOfPoints(50)
    cline.AddFieldFunction(fflist[0])
    cline.Initialize()
    cline.Execute()

    # Volume integrations
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)

    ffvol = FieldFunctionInterpolationVolume()
    ffvol.SetOperationType("max")
    ffvol.SetLogicalVolume(vol0)
    ffvol.AddFieldFunction(fflist[0])
    ffvol.Initialize()
    ffvol.Execute()
    maxval = ffvol.GetValue()
    if rank == 0:
        print(f"Max-value={maxval:.6f}")
