#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1-group infinite-medium transient with a ramped source and analytic solution.

ODE (v=1): d(phi)/dt + sigma_t * phi = Q(t)
Q(t) ramps linearly from 0 at t=0 to Q0 at t=t_ramp, then stays at Q0.
"""

import os
import sys
import math

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


def ramp_q(time_value: float, q0: float, t_ramp: float) -> float:
    if time_value <= 0.0:
        return 0.0
    if time_value < t_ramp:
        return q0 * time_value / t_ramp
    return q0


def analytic_phi(time_value: float,
                 q0: float,
                 t_ramp: float,
                 sigma_t: float,
                 v: float) -> float:
    lam = v * sigma_t
    if time_value <= t_ramp:
        a = q0 / t_ramp
        return v * a * (
            time_value / lam
            - 1.0 / (lam * lam)
            + math.exp(-lam * time_value) / (lam * lam)
        )
    # value at t_ramp
    a = q0 / t_ramp
    phi_tr = v * a * (
        t_ramp / lam - 1.0 / (lam * lam) + math.exp(-lam * t_ramp) / (lam * lam)
    )
    dt = time_value - t_ramp
    return phi_tr * math.exp(-lam * dt) + (v * q0 / lam) * (1.0 - math.exp(-lam * dt))


if __name__ == "__main__":
    meshgen = FromFileMeshGenerator(filename="../../../../assets/mesh/cube3.2.msh")
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    sigma_t = 1.0
    v = 1.0
    q0 = 1.0
    t_ramp = 0.5

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, 0.0, v)

    def source_func(group: int, time_value: float) -> float:
        return ramp_q(time_value, q0, t_ramp)
    vol_src = VolumetricSource(block_ids=[0], strength_function=source_func)

    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=16, scattering_order=0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[vol_src],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={"verbose_inner_iterations": False},
    )

    solver = TransientSolver(problem=phys, initial_state="zero")
    solver.Initialize()
    solver.SetTheta(0.5)

    dt = 0.01
    stop_time = 1.0
    current_time = 0.0

    while current_time < stop_time:
        target_time = min(current_time + dt, stop_time)
        solver.SetTimeStep(target_time - current_time)
        solver.Advance()
        current_time = target_time

    fflist = phys.GetScalarFluxFieldFunction()
    monitor_volume = RPPLogicalVolume(infx=True, infy=True, infz=True)
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType("max")
    field_interp.SetLogicalVolume(monitor_volume)
    field_interp.AddFieldFunction(fflist[0])
    field_interp.Initialize()
    field_interp.Execute()
    phi_num = field_interp.GetValue()

    phi_exact = analytic_phi(stop_time, q0, t_ramp, sigma_t, v)
    rel_err = abs(phi_num - phi_exact) / phi_exact
    pass_flag = 1 if rel_err < 0.01 else 0

    if rank == 0:
        print(f"RAMP_SOURCE_ANALYTIC phi_num {phi_num:.6f} phi_exact {phi_exact:.6f}")
        print(f"RAMP_SOURCE_ANALYTIC_PASS {pass_flag}")
