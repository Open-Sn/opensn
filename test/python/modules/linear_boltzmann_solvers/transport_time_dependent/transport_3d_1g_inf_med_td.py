#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, TimeDependentSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

#    nodes = []
#    N = 10
#    L = 3.2
#    xmin = -L / 2
#    dx = L / N
#    for i in range(N + 1):
#        nodes.append(xmin + i * dx)
#    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    meshgen = FromFileMeshGenerator(filename="../../../../assets/mesh/cube3.2.msh")
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    num_groups = 1 
    xs_diag = MultiGroupXS()
    xs_diag.LoadFromOpenSn("simple_1g_td.cxs")

    strength = [0.0 for _ in range(num_groups)]
    strength[0] = 122.58/(3.2*3.2*3.2)
    mg_src = VolumetricSource(block_ids=[0],
                              group_strength=strength,
                              start_time=0.0,
                              end_time=1.0)

    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=16, scattering_order=0)

    gs0 = [0, num_groups - 1]
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": gs0,
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 500,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_diag},
        ],
        volumetric_sources=[mg_src],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        time_dependent=True,
    )

    solver = TimeDependentSourceSolver(problem=phys, dt=0.05, theta=1.0, stop_time=1.0)
    solver.Initialize()
    solver.Execute()

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
        print(f"Peak flux = {flux_max:.6f}")
