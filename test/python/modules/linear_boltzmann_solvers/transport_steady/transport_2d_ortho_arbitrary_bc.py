#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D orthogonal mesh, 1g, pure absorber. Arbirary boundary condition on xmin, applied only to
directions with positive omega.x and omega.y. All other boundaries are vacuum.
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
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    num_procs = 4

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    nodes = []
    N = 64
    L = 1.0
    xmin = 0.0
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_mat = MultiGroupXS()
    xs_mat.CreateSimpleOneGroup(1.0, 0.0)
    num_groups = xs_mat.num_groups

    pquad = GLCProductQuadrature2DXY(n_polar=16, n_azimuthal=32, scattering_order=0)

    # Mask: angles where omega.x > 0 and omega.y > 0
    num_angles = len(pquad.omegas)
    angle_mask = np.zeros(num_angles, dtype=bool)
    norm = 0.0
    sum = 0.0

    for i, (omega, w) in enumerate(zip(pquad.omegas, pquad.weights)):
        ox = omega.x
        oy = omega.y
        oz = omega.z
        if ox > 0.0 and oy > 0.0:
            angle_mask[i] = True
            sum += ox * w * num_groups

    if sum <= 0.0:
        raise RuntimeError("xmin_bc_func: no angles satisfy ox > 0 and oy > 0")

    norm_factor = sum * L

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
        ],
    )

    # Initialize and Execute Solver
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    fflist = phys.GetScalarFieldFunctionList()
    monitor_volume = RPPLogicalVolume(infx=True, infy=True, infz=True)
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType("max")
    field_interp.SetLogicalVolume(monitor_volume)
    field_interp.AddFieldFunction(fflist[0])
    field_interp.Initialize()
    field_interp.Execute()
    flux_max = field_interp.GetValue()

    if rank == 0:
        print(f"Max phi = {flux_max:.6f}")
