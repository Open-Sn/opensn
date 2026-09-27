#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Exact theta-scheme recurrence for a fissile infinite medium with precursors.

A 1D slab with reflecting boundaries on both ends (an infinite medium) holds a
subcritical 1-group material with scattering, prompt fission, and two
delayed-neutron precursor families, driven by the ramp source Q(t) = 1 + 2t
from a zero initial state. For each theta in {1, 1/2, 0.7} and both sweep types,
the transport solution must match the scalar theta-scheme recurrence (see the
theory manual, Time Discretization), evaluated independently below:

  (phi^{n+th} - phi^n)/(v th dt) + sigma_a phi^{n+th}
      = nu_p sigma_f phi^{n+th} + sum_j lambda_j C_j^{n+th} + Q(t^n + th dt),
  C_j^{n+th} = (C_j^n + th dt gamma_j nu_d sigma_f phi^{n+th}) / (1 + th dt lambda_j),
  x^{n+1} = (x^{n+th} - (1 - th) x^n) / th   for x = phi, C_j.

This checks that fission within the groupset is implicit, that the precursor
update is consistent with the delayed source for any theta, and that sources
are evaluated at t^{n+theta}. It also checks the exact discrete balance
(N^{n+1} - N^n)/dt = rates at t^{n+theta} through the balance table.
Finally, a source callback that fails at t^{n+theta} verifies that a failed
advance restores the problem time to t^n.

Checks:
  THETA_EXACT_MAX_REL_ERR: max relative error in phi over all steps and cases.
  THETA_BALANCE_MAX_REL_RESID: max |inventory residual| / |inventory change|.
Both are at the level of the inner solver tolerance (1e-12).
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
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def source_strength(group, time):
    return 1.0 + 2.0 * time


def failing_source_strength(group, time):
    if time > 0.0:
        raise RuntimeError("intentional source failure")
    return 1.0


def reference(theta, dt, num_steps):
    # xs1g_delayed_deep_sub_2p_decay_test.cxs
    v, sigma_a = 1.0, 1.0 - 0.7
    nu_p_sigma_f, nu_d_sigma_f = 1.987 * 0.05, 0.013 * 0.05
    lam, gamma = [0.1, 1.0], [0.6, 0.4]

    phi, conc, history = 0.0, [0.0, 0.0], []
    t = 0.0
    for _ in range(num_steps):
        h = theta * dt
        # Implicit delayed coefficient and old-inventory decay term
        a_d = sum(lam[j] * h * gamma[j] * nu_d_sigma_f / (1.0 + h * lam[j]) for j in range(2))
        s_d = sum(lam[j] * conc[j] / (1.0 + h * lam[j]) for j in range(2))
        tau = 1.0 / (v * h)
        rhs = tau * phi + s_d + source_strength(0, t + h)
        phi_th = rhs / (tau + sigma_a - nu_p_sigma_f - a_d)
        conc_th = [(conc[j] + h * gamma[j] * nu_d_sigma_f * phi_th) / (1.0 + h * lam[j])
                   for j in range(2)]
        phi = (phi_th - (1.0 - theta) * phi) / theta
        conc = [(conc_th[j] - (1.0 - theta) * conc[j]) / theta for j in range(2)]
        t += dt
        history.append(phi)
    return history


def max_phi(problem):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("max")
    ffi.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    ffi.AddFieldFunction(problem.GetScalarFluxFieldFunction()[0])
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    nodes = [i * 0.25 for i in range(5)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "../../../../assets/xs/xs1g_delayed_deep_sub_2p_decay_test.cxs"))

    dt, num_steps = 0.2, 8
    max_err, max_resid = 0.0, 0.0
    for sweep_type in ("AAH", "CBC"):
        for theta in (1.0, 0.5, 0.7):
            problem = DiscreteOrdinatesProblem(
                mesh=grid,
                num_groups=1,
                time_dependent=True,
                sweep_type=sweep_type,
                groupsets=[
                    {
                        "groups_from_to": (0, 0),
                        "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4,
                                                                        scattering_order=0),
                        "inner_linear_method": "petsc_gmres",
                        "l_abs_tol": 1.0e-12,
                        "l_max_its": 500,
                    }
                ],
                xs_map=[{"block_ids": [0], "xs": xs}],
                volumetric_sources=[VolumetricSource(block_ids=[0],
                                                     strength_function=source_strength)],
                boundary_conditions=[
                    {"name": "zmin", "type": "reflecting"},
                    {"name": "zmax", "type": "reflecting"},
                ],
                options={
                    "save_angular_flux": True,
                    "use_precursors": True,
                    "verbose_inner_iterations": False,
                },
            )
            solver = TransientSolver(problem=problem, dt=dt, stop_time=1.0e9, theta=theta,
                                     initial_state="zero", verbose=False)
            solver.Initialize()
            ref = reference(theta, dt, num_steps)
            for n in range(num_steps):
                solver.Advance()
                max_err = max(max_err, abs(max_phi(problem) / ref[n] - 1.0))
                table = solver.ComputeBalanceTable()
                resid = table["inventory_residual"] / table["actual_inventory_change"]
                max_resid = max(max_resid, abs(resid))

                # The balance belongs to the completed step even if the next step changes theta.
                solver.SetTheta(0.5 if theta == 1.0 else 1.0)
                table = solver.ComputeBalanceTable()
                solver.SetTheta(theta)
                resid = table["inventory_residual"] / table["actual_inventory_change"]
                max_resid = max(max_resid, abs(resid))

    failure_problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        time_dependent=True,
        sweep_type="AAH",
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0],
                                             strength_function=failing_source_strength)],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={"save_angular_flux": True, "verbose_inner_iterations": False},
    )
    failure_solver = TransientSolver(problem=failure_problem, dt=dt, stop_time=dt, theta=0.5,
                                     initial_state="zero", verbose=False)
    failure_solver.Initialize()
    try:
        failure_solver.Advance()
        raise RuntimeError("advance unexpectedly succeeded with a failing source")
    except RuntimeError as error:
        if "intentional source failure" not in str(error):
            raise
    if failure_problem.GetTime() != 0.0:
        raise RuntimeError("failed advance did not restore the problem time")

    if rank == 0:
        print(f"THETA_EXACT_MAX_REL_ERR {max_err:.6e}")
        print(f"THETA_BALANCE_MAX_REL_RESID {max_resid:.6e}")
