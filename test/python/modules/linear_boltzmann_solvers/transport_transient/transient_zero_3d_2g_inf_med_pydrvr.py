#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# 2-group, infinite-medium transient with downscatter (group 0 -> group 1).
# 3.2 cm reflecting cube (infinite medium) with:
# g0 (fast):    sigma_t0 = 1.0 cm^-1, v0 = 2.0 cm/s
# g1 (thermal): sigma_t1 = 0.8 cm^-1, v1 = 0.5 cm/s
# sigma_s(0 -> 1) = 0.5 cm^-1, all other sigma_s = 0
# Constant in time source in g0 only:
# Q0 = 122.58 / 3.2^3 cm^-3 s^-1, for all t >= 0
#
# (1/v0) d(phi0)/dt + sigma_t0 * phi0 = Q0
# (1/v1) d(phi1)/dt + sigma_t1 * phi1 = sigma_s(0->1) * phi0(t)
# phi0(t) = (Q0 / sigma_t0) * (1 - exp(-v0 * sigma_t0 * t)) = Q0 * (1 - exp(-2 t))
# phi1(t) = exp(-v1 * sigma_t1 * t) *
#           [ v1 * sigma_s(0->1) * integral_0^t exp(v1 * sigma_t1 * s) * phi0(s) ds ]
# phi0(1s) ~= 3.235,  phi1(1s) ~= 0.458
# phi0(2s) ~= 3.672,  phi1(2s) ~= 1.036

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
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    meshgen = FromFileMeshGenerator(filename="../../../../assets/mesh/cube3.2.msh")
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    xs_diag = MultiGroupXS()
    xs_diag.LoadFromOpenSn(
        os.path.join(os.path.dirname(__file__), "simple_2g_downscatter_td.cxs")
    )
    num_groups = xs_diag.num_groups

    # Total source in group 0, converted to volumetric rate
    Q_tot = 122.58
    Q_vol = Q_tot / (3.2 * 3.2 * 3.2)

    strength = [0.0 for _ in range(num_groups)]
    strength[0] = Q_vol  # source only in group 0
    strength[1] = 0.0

    # Volumetric source is effectively always on
    mg_src = VolumetricSource(block_ids=[0], group_strength=strength)

    # Angular quadrature
    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=16, scattering_order=0)

    gs0 = [0, num_groups - 1]

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        time_dependent=True,
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
        options={
            "verbose_inner_iterations": False,
        },
    )

    # Create the time-dependent solver without stop_time, we will loop in Python
    solver = TransientSolver(problem=phys, verbose=False, initial_state="zero")
    solver.Initialize()

    # Time stepping parameters (constant dt)
    dt = 0.05
    theta = 0.5
    stop_time = 2.0
    current_time = 0.0
    step = 0
    solver.SetTheta(theta)

    fflist = phys.GetScalarFieldFunctionList()
    monitor_volume = RPPLogicalVolume(infx=True, infy=True, infz=True)

    while current_time < stop_time:
        target_time = min(current_time + dt, stop_time)
        step_dt = target_time - current_time

        # dt is constant here with the exception of the last step.
        # We adjust dt for the last step so that we get the solution
        # exactly at stop_time
        solver.SetTimeStep(step_dt)

        if rank == 0:
            print("")
            print(
                f"*************** Time step #{step:d}  t = {target_time:.6f} "
                f"(from {current_time:.6f}, dt = {step_dt:.6f}, theta = {theta:.3f}) "
                f"***************"
            )

        # Advance the solution from current_time to target_time
        solver.Advance()

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
            print("Max phi0 = {:.6f}".format(flux_max_g0))
            print("Max phi1 = {:.6f}".format(flux_max_g1))

        current_time = target_time
        step += 1
