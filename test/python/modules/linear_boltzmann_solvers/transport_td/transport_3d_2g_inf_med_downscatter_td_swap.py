#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# 2-group, infinite-medium transient with downscatter (group 0 -> group 1).
# Cross sections are swapped at t=0.5 s:
#   before: sigma_t0=1.0, sigma_t1=0.8, sigma_s01=0.5
#   after : sigma_t0=2.0, sigma_t1=1.2, sigma_s01=0.6
# Velocities remain v0=2.0 cm/s, v1=0.5 cm/s.
#
# (1/v0) d(phi0)/dt + sigma_t0 * phi0 = Q0
# (1/v1) d(phi1)/dt + sigma_t1 * phi1 = sigma_s01 * phi0(t)
# With backward Euler, dt=0.05: phi0(1s) ~= 1.939573, phi1(1s) ~= 0.384769

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

    xs_diag = MultiGroupXS()
    xs_diag.LoadFromOpenSn("simple_2g_downscatter_td.cxs")
    xs_diag_swap = MultiGroupXS()
    xs_diag_swap.LoadFromOpenSn("simple_2g_downscatter_td_swap.cxs")
    num_groups = xs_diag.num_groups

    Q_tot = 122.58
    Q_vol = Q_tot / (3.2 * 3.2 * 3.2)

    strength = [0.0 for _ in range(num_groups)]
    strength[0] = Q_vol   # source only in group 0
    strength[1] = 0.0

    mg_src = VolumetricSource(block_ids=[0],
                              group_strength=strength,
                              start_time=0.0,
                              end_time=1.0e9)  # effectively always on

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
                "inner_linear_method": "petsc_gmres",
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
    )

    solver = TimeDependentSourceSolver(problem=phys, verbose=False)
    solver.Initialize()

    dt = 0.05
    theta = 1.0
    step = 0
    stop_time = 1.0
    swap_time = 0.5
    current_time = 0.0
    swapped = False
    solver.SetTheta(theta)

    while current_time < stop_time:
        target_time = min(current_time + dt, stop_time)
        step_dt = target_time - current_time
        solver.SetTimeStep(step_dt)

        if rank == 0:
            print("")
            print(
                f"*************** Time step #{step:d}  t = {target_time:.6f} "
                f"(from {current_time:.6f}, dt = {step_dt:.6f}, theta = {theta:.3f}) "
                f"***************"
            )

        solver.Advance()

        if (not swapped) and target_time >= swap_time:
            phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_diag_swap}])
            swapped = True

        current_time = target_time
        step = step + 1

    fflist = phys.GetScalarFieldFunctionList()
    monitor_volume = RPPLogicalVolume(infx=True, infy=True, infz=True)

    # Group 0
    ff_interp_g0 = FieldFunctionInterpolationVolume()
    ff_interp_g0.SetOperationType("max")
    ff_interp_g0.SetLogicalVolume(monitor_volume)
    ff_interp_g0.AddFieldFunction(fflist[0])
    ff_interp_g0.Initialize()
    ff_interp_g0.Execute()
    flux_max_g0 = ff_interp_g0.GetValue()

    # Group 1
    ff_interp_g1 = FieldFunctionInterpolationVolume()
    ff_interp_g1.SetOperationType("max")
    ff_interp_g1.SetLogicalVolume(monitor_volume)
    ff_interp_g1.AddFieldFunction(fflist[1])
    ff_interp_g1.Initialize()
    ff_interp_g1.Execute()
    flux_max_g1 = ff_interp_g1.GetValue()

    if rank == 0:
        print("Max phi0(1s) = {:.6f}".format(flux_max_g0))
        print("Max phi1(1s) = {:.6f}".format(flux_max_g1))
