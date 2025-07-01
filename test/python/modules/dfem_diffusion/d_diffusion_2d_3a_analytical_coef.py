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
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    # Setup mesh
    nodes = []
    N = 100
    L = 2
    xmin = -L / 2
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    def D_coef(i, pt):
        return 3.0 + pt.x + pt.y

    def Q_ext(i, pt):
        return pt.x * pt.x

    def Sigma_a(i, pt):
        return pt.x * pt.y * pt.y

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

    # DFEM solver
    phys = DFEMDiffusionSolver(
        name="DFEMDiffusionSolver",
        mesh=grid,
        residual_tolerance=1e-8,
    )
    phys.SetOptions(boundary_conditions=[
        {"boundary": e_bndry, "type": "dirichlet", "coeffs": [0.0]},
        {"boundary": n_bndry, "type": "dirichlet", "coeffs": [0.0]},
        {"boundary": s_bndry, "type": "dirichlet", "coeffs": [0.0]},
        {"boundary": w_bndry, "type": "dirichlet", "coeffs": [0.0]},
    ])
    phys.SetDCoefFunction(D_coef)
    phys.SetQExtFunction(Q_ext)
    phys.SetSigmaAFunction(Sigma_a)
    phys.Initialize()
    phys.Execute()

    # Get field functions
    fflist = phys.GetFieldFunctions()

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
