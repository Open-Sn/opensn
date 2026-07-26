#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D analytic test of multi-family delayed-neutron precursor tracking.

This test is a generalization of transient_keigen_3d_delayed_analytic (which uses a
single precursor family and a hand-derived quadratic "inhour equation") to two precursor
families with well-separated decay constants. The reference solution here is obtained by
directly diagonalizing the exact governing linear ODE system.

The problem is a homogeneous, fully-reflecting 1D slab. A PowerIterationKEigenSolver
converges an initial critical, 2-precursor-family state, giving flux phi_0 and a known
initial precursor inventory C_j(0) = (raw_yield_j / sum(raw_yield)) / lambda_j *
nu_delayed * sigma_f * phi_0.

At t=0, the cross-section map is swapped to a deeply subcritical material with the same two
precursor families (so the precursor inventory is correctly carried forward across the swap,
unlike a swap to a zero-precursor material, which we intentionally discard history for).
Because the medium is homogeneous and fully reflecting, the transport equation for the
post-swap system reduces exactly to the coupled linear ODE system, for state vector
y = [phi, C_1, C_2]^T:

  dy/dt = M @ y,  M = [[v*(nu_prompt*sigma_f - sigma_a),  v*lambda_1,  v*lambda_2],
                        [beta_1*nu_delayed*sigma_f,        -lambda_1,   0        ],
                        [beta_2*nu_delayed*sigma_f,         0,         -lambda_2 ]]

with sigma_a = sigma_t - sigma_s0 (in-group scattering does not remove particles in a 1-group
problem) and beta_j the normalized precursor yields. This has the exact closed-form solution
y(t) = V @ diag(exp(eig*t)) @ V^-1 @ y(0), obtained via eigendecomposition of M.

If a future regression drops the decay of precursor inventory carried over from the previous
time step the computed flux will show a fast, incorrect initial transient instead of tracking
this reference -- the test fails loudly (order-of-magnitude, not a marginal tolerance miss)
because the analytic solution's late-time behavior is dominated by the slowest precursor
family's decay.

DECAY_ANALYTIC_PASS is 1 if the space-time flux matches this eigendecomposition-based
closed-form solution within rel_tol for every recorded time step, and the late-time logarithmic
decay rate matches the dominant eigenvalue of M.
"""

import math
import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        TransientSolver,
    )
    from pyopensn.xs import MultiGroupXS


def read_value(path, block_name, index=0):
    return read_values(path, block_name)[index]


def read_values(path, block_name):
    begin = f"{block_name}_BEGIN"
    end = f"{block_name}_END"
    in_block = False
    values = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line == begin:
                in_block = True
                continue
            if line == end:
                break
            if in_block:
                parts = line.split()
                if len(parts) >= 2:
                    values.append(float(parts[1]))
    if not values:
        raise RuntimeError(f"Failed to find values for {block_name} in {path}")
    return values


def max_phi(problem):
    fflist = problem.GetScalarFluxFieldFunction()
    monitor_volume = RPPLogicalVolume(infx=True, infy=True, infz=True)
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType("max")
    field_interp.SetLogicalVolume(monitor_volume)
    field_interp.AddFieldFunction(fflist[0])
    field_interp.Execute()
    return field_interp.GetValue()


def build_rate_matrix(v, sigma_a, nu_prompt_sigma_f, betas, nu_delayed_sigma_f, decay_consts):
    """Builds the exact 0-D governing matrix M for dy/dt = M @ y, y = [phi, C_1, ..., C_J]."""
    n_precursors = len(decay_consts)
    dim = 1 + n_precursors
    matrix = np.zeros((dim, dim))
    matrix[0, 0] = v * (nu_prompt_sigma_f - sigma_a)
    for j in range(n_precursors):
        matrix[0, 1 + j] = v * decay_consts[j]
        matrix[1 + j, 0] = betas[j] * nu_delayed_sigma_f
        matrix[1 + j, 1 + j] = -decay_consts[j]
    return matrix


class ExactLinearResponse:
    """Exact closed-form solution of dy/dt = M @ y via eigendecomposition of M."""

    def __init__(self, matrix, y0):
        eigvals, eigvecs = np.linalg.eig(matrix)
        self.eigvals = eigvals
        self.eigvecs = eigvecs
        self.coeffs = np.linalg.solve(eigvecs, y0)
        # Dominant (largest real part) eigenvalue, for the independent late-time slope check.
        self.dominant_eigval = eigvals[np.argmax(eigvals.real)]

    def __call__(self, t):
        y = self.eigvecs @ (self.coeffs * np.exp(self.eigvals * t))
        imag_scale = max(abs(y.real).max(), 1.0e-300)
        if abs(y.imag).max() > 1.0e-6 * imag_scale:
            raise RuntimeError("Unexpectedly large imaginary component in analytic solution.")
        return y.real


if __name__ == "__main__":
    n_cells = 20
    length = 8.0
    dx = length / n_cells
    nodes = [i * dx for i in range(n_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    base_dir = os.path.join(os.path.dirname(__file__), "../../../../assets/xs")
    xs_crit_path = os.path.join(base_dir, "xs1g_delayed_crit_2p_decay_test.cxs")
    xs_sub_path = os.path.join(base_dir, "xs1g_delayed_deep_sub_2p_decay_test.cxs")

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(xs_crit_path)

    xs_sub = MultiGroupXS()
    xs_sub.LoadFromOpenSn(xs_sub_path)

    pquad = GLProductQuadrature1DSlab(n_polar=8, scattering_order=0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )

    keigen = PowerIterationKEigenSolver(problem=problem, max_iters=200, k_tol=1.0e-12)
    keigen.Initialize()
    keigen.Execute()

    phi_0 = max_phi(problem)
    if phi_0 <= 0.0:
        raise RuntimeError("Initial scalar flux must be positive.")

    # Precursor decay data and secular-equilibrium inventory at t=0, using the same formula
    # OpenSn's own ComputePrecursors uses.
    raw_yields = read_values(xs_crit_path, "PRECURSOR_FRACTIONAL_YIELDS")
    decay_consts = read_values(xs_crit_path, "PRECURSOR_DECAY_CONSTANTS")
    sum_raw_yields = sum(raw_yields)
    betas = [y / sum_raw_yields for y in raw_yields]
    nu_delayed_crit = read_value(xs_crit_path, "NU_DELAYED")
    sigma_f_crit = read_value(xs_crit_path, "SIGMA_F")

    precursor_c0 = [
        beta_j / lam_j * nu_delayed_crit * sigma_f_crit * phi_0
        for beta_j, lam_j in zip(betas, decay_consts)
    ]

    # Post-step (deeply subcritical) material properties for the governing matrix.
    v = read_value(xs_sub_path, "VELOCITY")
    sigma_t_sub = read_value(xs_sub_path, "SIGMA_T")
    # TRANSFER_MOMENTS entries are "M_GFROM_GTO_VAL gp g val"; the in-group (0 0 0) moment is the
    # 4th whitespace-separated token on that data line.
    sigma_s0_sub = None
    with open(xs_sub_path, "r", encoding="utf-8") as handle:
        in_block = False
        for line in handle:
            line = line.strip()
            if line == "TRANSFER_MOMENTS_BEGIN":
                in_block = True
                continue
            if line == "TRANSFER_MOMENTS_END":
                break
            if in_block and line.startswith("M_GFROM_GTO_VAL"):
                parts = line.split()
                sigma_s0_sub = float(parts[-1])
    if sigma_s0_sub is None:
        raise RuntimeError(f"Failed to find in-group transfer moment in {xs_sub_path}")
    sigma_a_sub = sigma_t_sub - sigma_s0_sub
    nu_prompt_sub = read_value(xs_sub_path, "NU_PROMPT")
    sigma_f_sub = read_value(xs_sub_path, "SIGMA_F")
    nu_delayed_sub = read_value(xs_sub_path, "NU_DELAYED")

    rate_matrix = build_rate_matrix(
        v=v,
        sigma_a=sigma_a_sub,
        nu_prompt_sigma_f=nu_prompt_sub * sigma_f_sub,
        betas=betas,
        nu_delayed_sigma_f=nu_delayed_sub * sigma_f_sub,
        decay_consts=decay_consts,
    )
    y0 = np.array([phi_0] + precursor_c0)
    reference = ExactLinearResponse(rate_matrix, y0)

    # Shut off fission entirely: NOT via a zero-precursor swap (OpenSn intentionally discards
    # precursor history when a cell's material has zero precursor families), but via a deeply
    # subcritical material that keeps the same two precursor families, so the inventory
    # computed above is correctly carried forward across the swap.
    problem.SetTimeDependentMode()
    solver = TransientSolver(problem=problem, verbose=False, initial_state="existing")
    solver.Initialize()
    problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_sub}])

    dt = 0.05
    n_steps = 1800
    rel_tol = 5.0e-3

    solver.SetTimeStep(dt)
    solver.SetTheta(0.5)

    max_rel_err = 0.0
    times = []
    phi_numeric = []
    phi_reference = []

    if rank == 0:
        print("step time phi_numeric phi_analytic")

    for step in range(1, n_steps + 1):
        solver.Advance()
        t_now = problem.GetTime()
        phi_num = max_phi(problem)
        phi_ref = reference(t_now)[0]

        times.append(t_now)
        phi_numeric.append(phi_num)
        phi_reference.append(phi_ref)

        rel_err = abs(phi_num - phi_ref) / max(abs(phi_ref), 1.0e-300)
        max_rel_err = max(max_rel_err, rel_err)

        if rank == 0:
            print(f"{step:4d} {t_now:10.4e} {phi_num:14.8e} {phi_ref:14.8e}")

    pointwise_ok = max_rel_err < rel_tol

    # Independent cross-check: at late times the response is dominated by the eigenmode of M
    # with the largest real part, so a log-linear fit over the final 10% of the run should
    # recover that eigenvalue. This uses a different mathematical method (linear regression on
    # log-flux) than the primary pointwise eigendecomposition comparison above.
    tail_start = 9 * n_steps // 10
    tail_t = times[tail_start:]
    tail_log_phi = [math.log(p) for p in phi_numeric[tail_start:]]
    n_tail = len(tail_t)
    mean_t = sum(tail_t) / n_tail
    mean_log_phi = sum(tail_log_phi) / n_tail
    cov = sum((tt - mean_t) * (lp - mean_log_phi) for tt, lp in zip(tail_t, tail_log_phi))
    var = sum((tt - mean_t) ** 2 for tt in tail_t)
    fitted_slope = cov / var
    expected_slope = reference.dominant_eigval.real
    slope_rel_err = abs(fitted_slope - expected_slope) / abs(expected_slope)
    slope_ok = slope_rel_err < 2.0e-2

    ok = pointwise_ok and slope_ok

    if rank == 0:
        print(f"MAX_REL_ERR {max_rel_err:.10e}")
        print(f"FITTED_LATE_TIME_SLOPE {fitted_slope:.10e}")
        print(f"EXPECTED_LATE_TIME_SLOPE {expected_slope:.10e}")
        print(f"DECAY_ANALYTIC_PASS {1 if ok else 0}")
