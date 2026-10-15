#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Solvers revalidate and reconfigure when the problem changes after Initialize().

All cases are a 1D slab with reflecting boundaries (an infinite medium), whose
solutions are exact:

1. Power iteration, then SetXSMap with the same cross sections, then Execute.
   The cross-section swap rebuilds the problem's WGS solvers; power iteration
   must reapply its source scopes (fission is its source, not part of the WGS
   operator). k must equal k_inf = nu sigma_f / (sigma_t - sigma_s) = 1.2.
2. Power iteration, then SetTimeDependentMode, then Execute must be rejected.
   After SetSteadyStateMode, Execute must give k_inf = 1.0.
3. Transient solver (pure absorber sigma_t = 1.5, v = 1, Q = 1, zero initial
   state, backward Euler), then SetSteadyStateMode, then Advance must be
   rejected. After SetTimeDependentMode, Advance must give the exact update
   phi_1 = v dt Q / (1 + v sigma_t dt).

Checks:
  PI_SWAP_K_ERR, PI_MODE_K_ERR: |k - k_inf|.
  PI_TIME_DEPENDENT_REJECTED, TRANSIENT_STEADY_REJECTED: 1 if rejected.
  TRANSIENT_STEP_MAX_REL_ERR: max |phi / phi_exact - 1| after the step.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        TransientSolver,
    )
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def make_problem(grid, xs, sources=None, time_dependent=False):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=sources or [],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
        time_dependent=time_dependent,
    )


def rejected(action, expected):
    try:
        action()
    except ValueError as error:
        return expected in str(error)
    return False


if __name__ == "__main__":
    xs_dir = os.path.join(os.path.dirname(__file__), "../../../../assets/xs")

    nodes = [i * 10.0 / 20 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    # Case 1: nu sigma_f = 2.0 * 0.18, sigma_a = 1.0 - 0.7.
    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(xs_dir, "xs1g_prompt_super.cxs"))
    problem = make_problem(grid, xs_super)
    pi = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-12, max_iters=500)
    pi.Initialize()
    problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])
    pi.Execute()
    pi_swap_err = abs(pi.GetEigenvalue() - 2.0 * 0.18 / 0.3)

    # Case 2: nu sigma_f = (1.987 + 0.013) * 0.15, sigma_a = 1.0 - 0.7. The cross sections have
    # velocities, so time-dependent mode is valid for the problem but not for power iteration.
    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(xs_dir, "xs1g_delayed_crit_2p.cxs"))
    problem = make_problem(grid, xs_crit)
    pi = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-12, max_iters=500)
    pi.Initialize()
    problem.SetTimeDependentMode()
    pi_td_rejected = rejected(pi.Execute, "Problem is in time-dependent mode")
    problem.SetSteadyStateMode()
    pi.Execute()
    pi_mode_err = abs(pi.GetEigenvalue() - 2.0 * 0.15 / 0.3)

    # Case 3
    sigma_t, velocity, source, dt = 1.5, 1.0, 1.0, 0.1
    xs_absorber = MultiGroupXS()
    xs_absorber.CreateSimpleOneGroup(sigma_t, 0.0, velocity)
    problem = make_problem(
        grid,
        xs_absorber,
        sources=[VolumetricSource(block_ids=[0], group_strength=[source])],
        time_dependent=True,
    )
    transient = TransientSolver(problem=problem, dt=dt, stop_time=1.0, theta=1.0,
                                initial_state="zero", verbose=False)
    transient.Initialize()
    problem.SetSteadyStateMode()
    transient_rejected = rejected(transient.Advance, "Problem is in steady-state mode")
    problem.SetTimeDependentMode()
    transient.Advance()
    phi_exact = velocity * dt * source / (1.0 + velocity * sigma_t * dt)
    phi = problem.GetPhiNewLocal()
    local_err = max(abs(value / phi_exact - 1.0) for value in phi) if len(phi) else 0.0
    transient_err = MPIAllReduce(local_err, "max")

    if rank == 0:
        print(f"PI_SWAP_K_ERR {pi_swap_err:.6e}")
        print(f"PI_TIME_DEPENDENT_REJECTED {int(pi_td_rejected)}")
        print(f"PI_MODE_K_ERR {pi_mode_err:.6e}")
        print(f"TRANSIENT_STEADY_REJECTED {int(transient_rejected)}")
        print(f"TRANSIENT_STEP_MAX_REL_ERR {transient_err:.6e}")
