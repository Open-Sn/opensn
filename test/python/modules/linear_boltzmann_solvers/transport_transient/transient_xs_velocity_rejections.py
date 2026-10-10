#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Time-dependent mode requires inverse velocities in every cross section.

A 1D slab with reflecting boundaries (an infinite medium) holds a pure absorber
(sigma_t = 1.5, v = 1) with a uniform source Q = 1, starting from zero. Backward
Euler (theta = 1) gives the exact infinite-medium update
  phi_{n+1} = (phi_n + v dt Q) / (1 + v sigma_t dt).

After the first step, replacing the cross sections with ones that have no
velocity data must be rejected, and the second step must still follow the
original material. Creating a time-dependent problem, or switching a
steady-state problem to time-dependent mode, with such cross sections must also
be rejected.

Checks:
  XS_VELOCITY_SWAP_REJECTED, XS_VELOCITY_CTOR_REJECTED,
  XS_VELOCITY_MODE_SWITCH_REJECTED: 1 if rejected with the velocity error.
  XS_VELOCITY_MAX_REL_ERR: max over both steps of |phi / phi_exact - 1|.
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
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

VELOCITY_ERROR = "requires VELOCITY or INV_VELOCITY data"


def make_problem(grid, xs, time_dependent):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 500,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={"save_angular_flux": True, "verbose_inner_iterations": False},
        time_dependent=time_dependent,
    )


def rejected(action):
    try:
        action()
    except ValueError as error:
        return VELOCITY_ERROR in str(error)
    return False


def max_rel_err(problem, phi_exact):
    phi = problem.GetPhiNewLocal()
    local = max(abs(value / phi_exact - 1.0) for value in phi) if len(phi) else 0.0
    return MPIAllReduce(local, "max")


if __name__ == "__main__":
    sigma_t, velocity, source, dt = 1.5, 1.0, 1.0, 0.1

    nodes = [i * 10.0 / 20 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, 0.0, velocity)
    xs_no_velocity = MultiGroupXS()
    xs_no_velocity.CreateSimpleOneGroup(3.0, 0.0)

    problem = make_problem(grid, xs, True)
    solver = TransientSolver(problem=problem, dt=dt, stop_time=1.0, theta=1.0,
                             initial_state="zero", verbose=False)
    solver.Initialize()

    phi_exact = 0.0
    max_err = 0.0
    swap_rejected = False
    for step in range(2):
        if step == 1:
            swap_rejected = rejected(
                lambda: problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_no_velocity}])
            )
        solver.Advance()
        phi_exact = (phi_exact + velocity * dt * source) / (1.0 + velocity * sigma_t * dt)
        max_err = max(max_err, max_rel_err(problem, phi_exact))

    ctor_rejected = rejected(lambda: make_problem(grid, xs_no_velocity, True))
    steady_problem = make_problem(grid, xs_no_velocity, False)
    mode_switch_rejected = rejected(steady_problem.SetTimeDependentMode)

    if rank == 0:
        print(f"XS_VELOCITY_SWAP_REJECTED {int(swap_rejected)}")
        print(f"XS_VELOCITY_CTOR_REJECTED {int(ctor_rejected)}")
        print(f"XS_VELOCITY_MODE_SWITCH_REJECTED {int(mode_switch_rejected)}")
        print(f"XS_VELOCITY_MAX_REL_ERR {max_err:.6e}")
