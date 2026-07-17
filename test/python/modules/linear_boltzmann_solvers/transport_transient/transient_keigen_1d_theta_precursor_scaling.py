#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D transient check: delayed fission source theta-sensitivity.

Compare delayed-vs-prompt FP ratio deltas at theta=1.0 and theta=0.5 for a
single step, starting from a critical, precursor-bearing steady state (so the
precursor inventory C_j(0) = beta_j/lambda_j * F_ss is nonzero at t=0).

The delayed-fission source has two contributions: (a) decay of the precursor
inventory carried in from the previous time step (~coeff * C_j(t_n), nearly
theta-insensitive since coeff = lambda_j/(1+eff_dt*lambda_j) barely changes
with eff_dt=theta*dt for eff_dt*lambda_j << 1), and (b) newly-born precursors
from fission during this step (~coeff * eff_dt * beta_j * F, which scales
close to linearly with eff_dt=theta*dt). For this test's parameters
(lambda=0.08, beta=0.0065, dt=1e-2), contribution (a) is roughly three orders
of magnitude larger than (b), so the theta-dependence of the newly-born term
is swamped and the ratio of deltas is expected to be very close to 1.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        TransientSolver,
    )
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator


def run_case(theta, use_precursors, xs):
    dx = 8.0 / 40
    nodes = [i * dx for i in range(40 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    pquad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": use_precursors,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )

    keigen = PowerIterationKEigenSolver(problem=phys)
    keigen.Initialize()
    keigen.Execute()

    phys.SetTimeDependentMode()

    solver = TransientSolver(problem=phys, initial_state="existing")
    solver.Initialize()
    solver.SetTheta(theta)
    solver.SetTimeStep(1.0e-2)

    fp0 = phys.ComputeFissionProduction("new")
    solver.Advance()
    fp1 = phys.ComputeFissionProduction("new")

    return fp1 / fp0


if __name__ == "__main__":
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(os.path.join(
        os.path.dirname(__file__),
        "../../../../assets/xs/xs1g_delayed_crit_1p.cxs")
    )

    ratio_delayed_t1 = run_case(1.0, True, xs)
    ratio_prompt_t1 = run_case(1.0, False, xs)
    ratio_delayed_t05 = run_case(0.5, True, xs)
    ratio_prompt_t05 = run_case(0.5, False, xs)

    delta_t1 = ratio_delayed_t1 - ratio_prompt_t1
    delta_t05 = ratio_delayed_t05 - ratio_prompt_t05

    ok = delta_t1 > 0.0 and delta_t05 > 0.0 and delta_t05 > 1.0e-10
    if ok:
        ratio = delta_t1 / delta_t05
    else:
        ratio = 0.0

    if rank == 0:
        print(f"DELAYED_THETA_RATIO {ratio:.12e}")
