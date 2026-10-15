#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
With `use_source_moments=True`, a solve requires the source moments to have been provided.

The source moments are the fixed source of such a problem and are provided after construction
(e.g. with ReadSourceMoments, which may also follow SteadyStateSourceSolver.Initialize(); see
transport_3d_1_poly_qmom_part2.py). Without them, a solve used to read an empty buffer. Every path
that reads them must instead raise a ValueError on every rank:

1. SteadyStateSourceSolver: Initialize() is accepted (the moments may still be loaded), Execute() is
   rejected.
2. TransientSolver on a time-dependent problem: Execute() (its first step) is rejected.
3. SetTimeDependentMode() on a steady-state problem without stored angular fluxes, which rebuilds
   them with the steady-state sources: rejected, and the problem stays in steady-state mode.

Checks: SOURCE_MOMENTS_REQUIRED_PASSED is the number of the 4 outcomes above (steady Initialize
accepted, steady Execute rejected, transient step rejected, mode switch rejected with the mode
unchanged) that hold on every rank.
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
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, TransientSolver
    from pyopensn.xs import MultiGroupXS

REQUIREMENT = "requires source moments"


def make_problem(grid, xs, time_dependent=False):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        time_dependent=time_dependent,
        options={
            "use_source_moments": True,
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )


def outcome(action):
    try:
        action()
    except ValueError as error:
        return "rejected" if REQUIREMENT in str(error) else "other"
    return "accepted"


if __name__ == "__main__":
    grid = OrthogonalMeshGenerator(node_sets=[[i * 0.1 for i in range(11)]]).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.5, 1.0)

    results = []

    # 1. Steady state: Initialize accepted, Execute rejected.
    steady = SteadyStateSourceSolver(problem=make_problem(grid, xs))
    results.append(outcome(steady.Initialize) == "accepted")
    results.append(outcome(steady.Execute) == "rejected")

    # 2. Transient step.
    transient = TransientSolver(problem=make_problem(grid, xs, time_dependent=True), dt=0.1,
                                theta=1.0, stop_time=0.1, initial_state="existing", verbose=False)
    transient.Initialize()
    results.append(outcome(transient.Execute) == "rejected")

    # 3. Switching to time-dependent mode rebuilds psi with the steady-state sources.
    switch_problem = make_problem(grid, xs)
    switch_rejected = outcome(switch_problem.SetTimeDependentMode) == "rejected"
    results.append(switch_rejected and not switch_problem.IsTimeDependent())

    passed = sum(int(MPIAllReduce(int(ok), "min")) for ok in results)
    if rank == 0:
        print(f"SOURCE_MOMENTS_REQUIRED_PASSED {passed}")
