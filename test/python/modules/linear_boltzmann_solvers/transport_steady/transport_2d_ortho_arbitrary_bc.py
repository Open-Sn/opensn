#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D orthogonal mesh on a rectangular domain with a void material.
Arbitrary boundary condition on xmin, applied only to directions with positive omega.x
and omega.y. All other boundaries are vacuum.
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

if __name__ == "__main__":

    num_procs = 4

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # 2D rectangular mesh
    Nx = 64
    Ny = 32
    Lx = 2.0
    Ly = 1.0
    xmin = 0.0
    ymin = 0.0
    dx = Lx / Nx
    dy = Ly / Ny
    nodes_x = []
    nodes_y = []
    for i in range(Nx + 1):
        nodes_x.append(xmin + i * dx)
    for i in range(Ny + 1):
        nodes_y.append(ymin + i * dy)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes_x, nodes_y])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Void
    xs_mat = MultiGroupXS()
    xs_mat.CreateSimpleOneGroup(0.0, 0.0)
    num_groups = xs_mat.num_groups

    # Angles with omega.x > 0 and omega.y > 0
    pquad = GLCProductQuadrature2DXY(n_polar=16, n_azimuthal=32, scattering_order=0)
    num_angles = len(pquad.omegas)
    angle_mask = np.zeros(num_angles, dtype=bool)
    sum = 0.0

    # Accumulate incoming current on xmin for all groups
    for i, (omega, w) in enumerate(zip(pquad.omegas, pquad.weights)):
        ox = omega.x
        oy = omega.y
        if ox > 0.0 and oy > 0.0:
            angle_mask[i] = True
            sum += ox * w * num_groups

    if sum <= 0.0:
        raise RuntimeError("xmin_bc_func: no angles satisfy ox > 0 and oy > 0")

    # Normalize so the total inflow through xmin integrates to 1.0
    norm_factor = sum * Ly

    # Return a constant inflow on xmin for the masked angles and all groups
    def xmin_bc_func(group_index, angle_index):
        if angle_mask[angle_index]:
            return 1.0 / norm_factor
        return 0.0

    xmin_bc = AngularFluxFunction(xmin_bc_func)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, num_groups - 1],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 1000,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_mat},
        ],
        boundary_conditions=[
            {"name": "xmin", "type": "arbitrary", "function": xmin_bc},
            {"name": "xmax", "type": "vacuum"},
            {"name": "ymin", "type": "vacuum"},
            {"name": "ymax", "type": "vacuum"},
        ],
        options={"save_angular_flux": True},
    )

    # Initialize and execute solver
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Compute leakage
    bndrys = ["xmin", "xmax", "ymin", "ymax"]
    leakage = phys.ComputeLeakage(bndrys)
    lkg_xmin = leakage["xmin"]
    lkg_xmax = leakage["xmax"]
    lkg_ymin = leakage["ymin"]
    lkg_ymax = leakage["ymax"]

    if rank == 0:
        print(f"Right leakage = {lkg_xmax.item()}")
        print(f"Left leakage = {lkg_xmin.item()}")
        print(f"Top leakage = {lkg_ymax.item()}")
        print(f"Bottom leakage = {lkg_ymin.item()}")
        print(f"Top + Right leakage = {(lkg_ymax + lkg_xmax).item()}")
