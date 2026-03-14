#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-step zero-init transient analytic ratio check versus density.

For dphi/dt + rho*sigma_t*phi = q with phi^0=0 and Backward Euler:
  phi^1 = q*dt / (1 + rho*sigma_t*dt).
Expected answers in this test:
  phi(rho=1) and phi(rho=2) each match the formula above, and
  phi2/phi1 = (1 + dt*sigma_t)/(1 + 2*dt*sigma_t) = 11/12 = 0.9166666667
for dt=0.1, sigma_t=1.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def max_phi(problem):
    ff = problem.GetScalarFluxFieldFunction()[0]
    lv = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("max")
    ffi.SetLogicalVolume(lv)
    ffi.AddFieldFunction(ff)
    ffi.Initialize()
    ffi.Execute()
    return ffi.GetValue()


def run_case(rho):
    nodes = [i / 20.0 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 1.0
    source_val = 1.5
    dt = 0.1

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, 0.0, 1.0)
    src = VolumetricSource(block_ids=[0], group_strength=[source_val], start_time=0.0, end_time=1.0)
    quad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        time_dependent=True,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": quad,
            "inner_linear_method": "petsc_richardson",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 200,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        density={"default_density": rho},
        options={"save_angular_flux": True},
    )

    solver = TransientSolver(problem=problem, dt=dt, theta=1.0, stop_time=dt, initial_state="zero")
    solver.Initialize()
    solver.Execute()

    phi_num = max_phi(problem)
    phi_analytic = source_val * dt / (1.0 + rho * sigma_t * dt)
    rel_err = abs(phi_num - phi_analytic) / phi_analytic

    return phi_num, rel_err


if __name__ == "__main__":
    phi1, err1 = run_case(1.0)
    phi2, err2 = run_case(2.0)
    ratio = phi2 / phi1 if phi1 != 0.0 else 0.0

    pass_flag = 1 if (err1 < 1.0e-6 and err2 < 1.0e-6 and abs(ratio - 0.9166666667) < 1.0e-6) else 0

    if rank == 0:
        print(f"DENSITY_TRANSIENT_ERR1={err1:.8e}")
        print(f"DENSITY_TRANSIENT_ERR2={err2:.8e}")
        print(f"DENSITY_TRANSIENT_RATIO={ratio:.8e}")
        print(f"DENSITY_TRANSIENT_PASS {pass_flag}")
