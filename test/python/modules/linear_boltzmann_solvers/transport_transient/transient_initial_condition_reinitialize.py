#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A steady-state initial condition (`read_initial_condition_path`) is loaded once per problem.

Calling TransientSolver.Initialize() again, or initializing a second transient solver on the same
problem, must succeed and continue from the problem's current state instead of reloading the file.

A 1 cm slab with reflecting boundaries on both ends is an infinite medium. The pure absorber
(sigma = 1 cm^-1, v = 1 cm/s) has a uniform source, so the flux is spatially constant and the
discrete solution is exact:

* The steady-state run with source Q1 writes phi0 = Q1 / sigma to the restart file.
* The transient problem has source Q2 and loads that file. Backward Euler then gives
      phi^{n+1} = (phi^n / (v dt) + Q2) / (sigma + 1 / (v dt)).
* Solver A is initialized twice, then runs N steps to t = N dt. Solver B (same problem, "existing"
  state) continues the problem's clock: its Initialize() leaves the time at N dt, and with the
  absolute stop_time (N + M) dt it runs M more steps. Reloading the initial condition would restart
  B from phi0, and resetting the clock to 0 would make it run N + M steps; both change the flux.
* Solver C ("zero" state) starts a new transient at t = 0.

Checks: REINITIALIZE_ACCEPTED and SECOND_SOLVER_ACCEPTED are 1; PHI_AFTER_A_REL_ERR and
PHI_AFTER_B_REL_ERR are max |phi - phi_exact| / phi_exact; RELOAD_WOULD_DIFFER and
CLOCK_RESET_WOULD_DIFFER are 1 if those defects would change the flux by more than 1e-3 (so the
test can detect them); B_INIT_TIME_ERR, B_FINAL_TIME_ERR and ZERO_START_TIME_ERR are the absolute
errors of the problem time after B.Initialize(), after B.Execute(), and after C.Initialize().
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, TransientSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def phi_extreme(phys, operation):
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType(operation)
    field_interp.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    field_interp.AddFieldFunction(phys.GetScalarFluxFieldFunction()[0])
    field_interp.Execute()
    return field_interp.GetValue()


def rel_err(phys, phi_exact):
    return max(abs(phi_extreme(phys, "max") - phi_exact),
               abs(phi_extreme(phys, "min") - phi_exact)) / phi_exact


def make_problem(grid, xs, q, options):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[q])],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            **options,
        },
    )


def backward_euler(phi, q, sigma, v, dt, steps):
    for _ in range(steps):
        phi = (phi / (v * dt) + q) / (sigma + 1.0 / (v * dt))
    return phi


def accepted(action):
    try:
        action()
    except ValueError:
        return False
    return True


if __name__ == "__main__":
    grid = OrthogonalMeshGenerator(node_sets=[[i * 0.1 for i in range(11)]]).Execute()
    grid.SetUniformBlockID(0)

    sigma, v = 1.0, 1.0
    q1, q2 = 2.0, 1.0
    dt, n_a, n_b = 0.1, 3, 2

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma, 0.0, v)

    restart_base = "transient_initial_condition_reinitialize"
    steady_problem = make_problem(
        grid, xs, q1, {"restart_writes_enabled": True, "write_restart_path": restart_base}
    )
    steady = SteadyStateSourceSolver(problem=steady_problem)
    steady.Initialize()
    steady.Execute()

    phys = make_problem(grid, xs, q2, {"read_initial_condition_path": restart_base})

    solver_a = TransientSolver(problem=phys, dt=dt, theta=1.0, stop_time=n_a * dt,
                               initial_state="existing", verbose=False)
    solver_a.Initialize()
    reinitialize_accepted = accepted(solver_a.Initialize)
    solver_a.Execute()
    phi0 = q1 / sigma
    phi_a = backward_euler(phi0, q2, sigma, v, dt, n_a)
    phi_after_a_err = rel_err(phys, phi_a)

    solver_b = TransientSolver(problem=phys, dt=dt, theta=1.0, stop_time=(n_a + n_b) * dt,
                               initial_state="existing", verbose=False)
    second_solver_accepted = accepted(solver_b.Initialize)
    b_init_time_err = abs(phys.GetTime() - n_a * dt)
    if second_solver_accepted:
        solver_b.Execute()
    b_final_time_err = abs(phys.GetTime() - (n_a + n_b) * dt)
    phi_b = backward_euler(phi_a, q2, sigma, v, dt, n_b)
    phi_after_b_err = rel_err(phys, phi_b)
    phi_b_reloaded = backward_euler(phi0, q2, sigma, v, dt, n_b)
    reload_would_differ = abs(phi_b_reloaded - phi_b) / phi_b > 1.0e-3
    phi_b_clock_reset = backward_euler(phi_a, q2, sigma, v, dt, n_a + n_b)
    clock_reset_would_differ = abs(phi_b_clock_reset - phi_b) / phi_b > 1.0e-3

    solver_c = TransientSolver(problem=phys, dt=dt, theta=1.0, stop_time=dt,
                               initial_state="zero", verbose=False)
    solver_c.Initialize()
    zero_start_time_err = abs(phys.GetTime())

    restart_file = f"{restart_base}{rank}.restart.h5"
    if os.path.exists(restart_file):
        os.remove(restart_file)

    if rank == 0:
        print(f"REINITIALIZE_ACCEPTED {int(reinitialize_accepted)}")
        print(f"SECOND_SOLVER_ACCEPTED {int(second_solver_accepted)}")
        print(f"PHI_AFTER_A_REL_ERR {phi_after_a_err:.6e}")
        print(f"PHI_AFTER_B_REL_ERR {phi_after_b_err:.6e}")
        print(f"RELOAD_WOULD_DIFFER {int(reload_would_differ)}")
        print(f"CLOCK_RESET_WOULD_DIFFER {int(clock_reset_would_differ)}")
        print(f"B_INIT_TIME_ERR {b_init_time_err:.6e}")
        print(f"B_FINAL_TIME_ERR {b_final_time_err:.6e}")
        print(f"ZERO_START_TIME_ERR {zero_start_time_err:.6e}")
