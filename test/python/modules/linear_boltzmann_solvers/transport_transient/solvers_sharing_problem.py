#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Solvers that share a problem do not affect each other's results.

All cases are a 1D slab with reflecting boundaries (an infinite medium), whose
solutions are exact. Material: sigma_t = 1, sigma_s = 0.7,
nu sigma_f = (1.987 + 0.013) * 0.1 = 0.2, so k_inf = 0.2 / 0.3 = 2/3 and the
fixed-source flux for Q = 1 is phi = Q / (0.3 - 0.2) = 10.

1. Fixed source, power iteration, fixed source, nonlinear k-eigenvalue, fixed
   source on one problem. The volumetric source is removed for the eigenvalue
   solves (an eigenvalue solver rejects a problem with external sources) and
   restored afterwards. Every fixed-source solve must give phi = 10 and every
   eigenvalue solve k = 2/3: an eigenvalue solver must not leave its source
   settings in the problem.
2. Two transient solvers with different time steps on one problem (pure
   absorber sigma_t = 1.5, v = 1, Q = 1, backward Euler), advanced
   alternately. Each step must give the exact update
   phi_{n+1} = (phi_n + v dt Q) / (1 + v sigma_t dt) with the stepping
   solver's own dt.

Checks:
  SHARED_EIGEN_SOURCE_REJECTED: 1 if power iteration rejects the sourced problem.
  SHARED_FIXED_SOURCE_MAX_REL_ERR: max |phi / 10 - 1| over the fixed-source solves.
  SHARED_K_MAX_ERR: max |k - 2/3| over the eigenvalue solves.
  SHARED_TRANSIENT_MAX_REL_ERR: max |phi / phi_exact - 1| over the transient steps.
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
        NonLinearKEigenSolver,
        PowerIterationKEigenSolver,
        SteadyStateSourceSolver,
        TransientSolver,
    )
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def make_problem(grid, xs, sources, time_dependent=False):
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
        volumetric_sources=sources,
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


def max_rel_err(problem, exact):
    phi = problem.GetPhiNewLocal()
    local = max(abs(value / exact - 1.0) for value in phi) if len(phi) else 0.0
    return MPIAllReduce(local, "max")


if __name__ == "__main__":
    nodes = [i * 10.0 / 20 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    # Case 1
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(
        os.path.join(os.path.dirname(__file__), "../../../../assets/xs/xs1g_delayed_sub_1p.cxs")
    )
    k_exact = 0.2 / 0.3
    phi_exact = 1.0 / (0.3 - 0.2)

    def source():
        return VolumetricSource(block_ids=[0], group_strength=[1.0])

    problem = make_problem(grid, xs, [source()])
    fixed_err = 0.0
    k_err = 0.0

    def fixed_source_solve():
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()
        return max_rel_err(problem, phi_exact)

    fixed_err = max(fixed_err, fixed_source_solve())

    # The documented false flag is a no-op. The source must still make the k-eigenvalue solver
    # reject this problem.
    problem.SetVolumetricSources(clear_volumetric_sources=False)
    pi = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-12, max_iters=1000)
    source_rejected = False
    try:
        pi.Initialize()
    except ValueError as error:
        source_rejected = "volumetric sources" in str(error)

    for eigen_solver in (
        PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-12, max_iters=1000),
        NonLinearKEigenSolver(problem=problem, nl_abs_tol=1.0e-12, nl_rel_tol=1.0e-12,
                              num_initial_power_iterations=5),
    ):
        problem.SetVolumetricSources(clear_volumetric_sources=True)
        eigen_solver.Initialize()
        eigen_solver.Execute()
        k_err = max(k_err, abs(eigen_solver.GetEigenvalue() - k_exact))
        problem.SetVolumetricSources(volumetric_sources=[source()])
        fixed_err = max(fixed_err, fixed_source_solve())

    # Case 2
    sigma_t, velocity, q = 1.5, 1.0, 1.0
    absorber = MultiGroupXS()
    absorber.CreateSimpleOneGroup(sigma_t, 0.0, velocity)
    problem = make_problem(
        grid,
        absorber,
        [VolumetricSource(block_ids=[0], group_strength=[q])],
        time_dependent=True,
    )
    solvers = [
        TransientSolver(problem=problem, dt=dt, stop_time=10.0, theta=1.0,
                        initial_state="zero", verbose=False)
        for dt in (0.1, 0.05)
    ]
    for solver in solvers:
        solver.Initialize()
    transient_err = 0.0
    phi = 0.0
    for step in range(4):
        dt = (0.1, 0.05)[step % 2]
        solvers[step % 2].Advance()
        phi = (phi + velocity * dt * q) / (1.0 + velocity * sigma_t * dt)
        transient_err = max(transient_err, max_rel_err(problem, phi))

    if rank == 0:
        print(f"SHARED_EIGEN_SOURCE_REJECTED {int(source_rejected)}")
        print(f"SHARED_FIXED_SOURCE_MAX_REL_ERR {fixed_err:.6e}")
        print(f"SHARED_K_MAX_ERR {k_err:.6e}")
        print(f"SHARED_TRANSIENT_MAX_REL_ERR {transient_err:.6e}")
