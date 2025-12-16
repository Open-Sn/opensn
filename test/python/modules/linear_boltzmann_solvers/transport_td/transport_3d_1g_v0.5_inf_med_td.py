#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 1-group, infinite-medium, transient in a 3.2 cm cube (all reflecting).
# The cube is a pure absorber with sigma_t = 1.0 cm^-1, sigma_s = 0, v = 0.5 cm/s
# and a total Q_tot = 122.58 particles/s on for t=[0, 1] s, then 0 for t>1 s.
# V = 3.2^3 cm^3, so volumetric Q = Q_tot / V ~= 3.7408 cm^-3 s^-1.
#
#   For 0 <= t <= 1: phi(t) = Q * (1 - e^{-t})
#   For t >= 1:      phi(t) = phi(1) * e^{-(t - 1)}
#   phi(1s) ~= 1.472
#   phi(2s) ~= 0.893

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

    meshgen = FromFileMeshGenerator(filename="../../../../assets/mesh/cube3.2.msh")
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    num_groups = 1
    xs_diag = MultiGroupXS()
    xs_diag.CreateSimpleOneGroup(1.0, 0.0, 0.5)

    strength = [0.0 for _ in range(num_groups)]
    strength[0] = 122.58 / (3.2 * 3.2 * 3.2)
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

    solver = TimeDependentSourceSolver(problem=phys, dt=0.05, theta=BackwardEuler, stop_time=2.0)
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
        print(f"Max phi(2s) = {flux_max:.6f}")
