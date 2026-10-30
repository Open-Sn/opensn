#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
IsConverged() of the steady-state, transient and nonlinear k-eigenvalue solvers, and
GetNumSweeps() of the nonlinear k-eigenvalue solver.

Medium: one group, sigma_t = 1, sigma_s = 0.9, on a 1D slab of length 10 cm with 20 cells,
vacuum faces and GL S8. A uniform source drives the fixed-source and transient problems; the
k-eigenvalue problem adds nu sigma_f = 0.15.

Each solver runs once with its default limits, which converge, and once with an iteration limit
that cannot be met:
  steady:    classic Richardson with l_max_its = 2 and l_abs_tol = 1e-12. Source iteration
             reduces the error by at most 0.9 per sweep, so two sweeps cannot reach 1e-12.
             PETSc Richardson with l_max_its = 1 performs a single sweep without a convergence
             test, which must not be reported as converged.
  transient: classic Richardson with l_max_its = 30 and backward Euler. A step of dt = 1e-4 adds
             the time absorption 1/(v dt) = 1e4, so the scattering ratio is about 9e-5 and the
             step converges in a few iterations. A step of dt = 1e3 leaves it at about 0.9,
             and 0.9^30 = 0.04 is far from 1e-12. The solver takes two small steps, then a
             large one, then a small one; Initialize() is then called again before one more
             small step.
  NLKE:      nl_max_its = 1.

Checks (1 = expected outcome, 0 = not):
  STEADY_CONVERGED, STEADY_LIMIT_NOT_CONVERGED, STEADY_NOT_CONVERGED_BEFORE_EXECUTE,
  STEADY_SINGLE_SWEEP_NOT_CONVERGED
  TRANSIENT_CONVERGED (after the two small steps), TRANSIENT_LIMIT_NOT_CONVERGED (after the
  large step), TRANSIENT_FAILURE_PERSISTS (after the next small step),
  TRANSIENT_NOT_CONVERGED_BEFORE_STEP (after the second Initialize()),
  TRANSIENT_RESET_BY_INITIALIZE (after the last small step)
  NLKE_CONVERGED, NLKE_LIMIT_NOT_CONVERGED
  NLKE_SWEEPS_COUNTED: GetNumSweeps() is positive and grows with a second Execute().
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
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.solver import TransientSolver, NonLinearKEigenSolver


def load_xs(fission):
    path = f"solver_is_converged_{int(fission)}_r{rank}.xs"
    with open(path, "w") as f:
        f.write("NUM_GROUPS 1\nNUM_MOMENTS 1\n"
                "SIGMA_T_BEGIN\n0 1.0\nSIGMA_T_END\n"
                "INV_VELOCITY_BEGIN\n0 1.0\nINV_VELOCITY_END\n"
                "TRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 0.9\nTRANSFER_MOMENTS_END\n")
        if fission:
            f.write("SIGMA_F_BEGIN\n0 0.06\nSIGMA_F_END\nNU_BEGIN\n0 2.5\nNU_END\n"
                    "CHI_BEGIN\n0 1.0\nCHI_END\n")
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    os.remove(path)
    return xs


def make_problem(grid, fission, l_max_its, method="classic_richardson"):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
                "inner_linear_method": method,
                "l_abs_tol": 1.0e-12,
                "l_max_its": l_max_its,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": load_xs(fission)}],
        volumetric_sources=[] if fission else [
            VolumetricSource(block_ids=[0], group_strength=[1.0])
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )


def steady(grid, l_max_its, execute=True, method="classic_richardson"):
    solver = SteadyStateSourceSolver(problem=make_problem(grid, False, l_max_its, method))
    solver.Initialize()
    if execute:
        solver.Execute()
    return solver.IsConverged()


if __name__ == "__main__":
    nodes = [i * 10.0 / 20 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    results = {
        "STEADY_CONVERGED": steady(grid, 1000),
        "STEADY_LIMIT_NOT_CONVERGED": not steady(grid, 2),
        "STEADY_NOT_CONVERGED_BEFORE_EXECUTE": not steady(grid, 1000, execute=False),
        "STEADY_SINGLE_SWEEP_NOT_CONVERGED": not steady(grid, 1, method="petsc_richardson"),
    }

    problem = make_problem(grid, False, 30)
    problem.SetTimeDependentMode()
    solver = TransientSolver(problem=problem, dt=1.0e-4, stop_time=1.0e9, theta=1.0,
                             initial_state="zero", verbose=False)
    solver.Initialize()
    solver.Advance()
    solver.Advance()
    results["TRANSIENT_CONVERGED"] = solver.IsConverged()
    solver.SetTimeStep(1.0e3)
    solver.Advance()
    results["TRANSIENT_LIMIT_NOT_CONVERGED"] = not solver.IsConverged()
    solver.SetTimeStep(1.0e-4)
    solver.Advance()
    results["TRANSIENT_FAILURE_PERSISTS"] = not solver.IsConverged()
    solver.Initialize()
    results["TRANSIENT_NOT_CONVERGED_BEFORE_STEP"] = not solver.IsConverged()
    solver.Advance()
    results["TRANSIENT_RESET_BY_INITIALIZE"] = solver.IsConverged()

    keigen_problem = make_problem(grid, True, 1000)
    nlke = NonLinearKEigenSolver(problem=keigen_problem, nl_abs_tol=1.0e-10, nl_max_its=50)
    nlke.Initialize()
    nlke.Execute()
    results["NLKE_CONVERGED"] = nlke.IsConverged()
    sweeps = nlke.GetNumSweeps()
    nlke.Execute()
    results["NLKE_SWEEPS_COUNTED"] = sweeps > 0 and nlke.GetNumSweeps() > sweeps

    nlke_limited = NonLinearKEigenSolver(problem=make_problem(grid, True, 1000),
                                         nl_abs_tol=1.0e-14, nl_max_its=1)
    nlke_limited.Initialize()
    nlke_limited.Execute()
    results["NLKE_LIMIT_NOT_CONVERGED"] = not nlke_limited.IsConverged()

    if rank == 0:
        for key, value in results.items():
            print(f"{key}={int(bool(value))}")
