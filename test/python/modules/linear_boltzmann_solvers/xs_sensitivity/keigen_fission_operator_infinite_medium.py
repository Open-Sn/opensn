#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
k-eigenvalue sensitivities use the fission production operator that the solvers apply.

A two-group 1D slab with reflecting boundaries is an infinite medium: the flux is flat and
isotropic, so the discrete problem reduces to k L v = F v with L = diag(sigma_t) - S0. With w the
left eigenvector (F^T w = k L^T w), the exact first-order sensitivities of k are

    dk = w^T (dF - k dL) v / (w^T L v),

so, relative to the coefficient x and to k,
  production (nu_sigma_f of source group g', spectra held fixed, so column g' of F scales):
      w^T F[:, g'] v_g' / (w^T F v),
  sigma_t of group g:            -sigma_t_g w_g v_g / (w^T L v),
  scatter S0[g][g'] (g' -> g):    S0[g][g'] w_g v_g' / (w^T L v).

The reference is evaluated with numpy from the same data. The scenario selects the fission data:
  matrix:      An explicit PRODUCTION_MATRIX that is not chi (x) nu_sigma_f, with a CHI block.
  delayed_on:  Prompt and delayed data with use_precursors=True. F is prompt plus delayed.
  delayed_off: The same data with use_precursors=False. F is prompt only.
  zero_column: An explicit PRODUCTION_MATRIX whose column 1 is zero, so nu_sigma_f[1] = 0. The
    column does not define the spectrum, so the relative production sensitivity to group 1 is
    exactly 0 and the absolute one raises ValueError. The exact relative sensitivity of group 1
    is 0, so the production error covers group 0 only. Group 1 neither fissions nor scatters up,
    so its importance w_1 is 0 and the exact sigma_t sensitivity of group 1 and the 1 <- 0 scatter
    sensitivity are 0 as well. An exact value of 0 is compared absolutely, relative to the largest
    sigma_t sensitivity.
Before the fix the postprocessor used chi (x) nu_sigma_f: matrix gave errors of 1e-3 to 1.5e-2 and
the delayed scenarios raised. The solves converge to 1e-12, so the acceptance tolerance is 1e-8.

With use_gpus=True the sweeps run on the GPU.

XS_SENS_K_SCALING_GUARDED is 1 when ApplyKEigenvalueScaling raises RuntimeError before Execute and
when the values are already scaled, leaves them unchanged then, and scales them the same way again
after a new Execute.
"""

import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.post import CrossSectionSensitivityPostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier


SIGMA_T = np.array([1.0, 1.5])
# S0[g][gp] is the P0 transfer from group gp to group g.
S0 = np.array([[0.5, 0.0], [0.3, 1.2]])
SIGMA_F = np.array([0.1, 0.2])
# Explicit production matrix F[g][gp] (gp -> g), deliberately not separable.
F_MATRIX = np.array([[0.40, 0.90], [0.05, 0.10]])
# Explicit production matrix with no production from group 1.
F_ZERO_COLUMN = np.array([[0.40, 0.0], [0.05, 0.0]])
CHI_PROMPT = np.array([0.9, 0.1])
NU_PROMPT_SIGMA_F = np.array([0.30, 0.70])
CHI_DELAYED = np.array([0.6, 0.4])
NU_DELAYED_SIGMA_F = np.array([0.003, 0.007])


def fission_operator(scenario):
    if scenario == "matrix":
        return F_MATRIX
    if scenario == "zero_column":
        return F_ZERO_COLUMN
    F = np.outer(CHI_PROMPT, NU_PROMPT_SIGMA_F)
    if scenario == "delayed_on":
        F = F + np.outer(CHI_DELAYED, NU_DELAYED_SIGMA_F)
    return F


def exact_relative_sensitivities(scenario):
    F = fission_operator(scenario)
    L = np.diag(SIGMA_T) - S0
    values, vectors = np.linalg.eig(np.linalg.solve(L, F))
    i = np.argmax(values.real)
    k, v = values[i].real, vectors[:, i].real
    values, vectors = np.linalg.eig(np.linalg.solve(L.T, F.T))
    w = vectors[:, np.argmax(values.real)].real
    production = [w @ (F[:, gp] * v[gp]) / (w @ F @ v) for gp in range(2)]
    sigma_t = [-SIGMA_T[g] * w[g] * v[g] / (w @ L @ v) for g in range(2)]
    scatter_1_from_0 = S0[1][0] * w[1] * v[0] / (w @ L @ v)
    return k, production, sigma_t, scatter_1_from_0


def write_xs(path, scenario):
    lines = ["NUM_GROUPS 2", "NUM_MOMENTS 1"]
    lines += ["SIGMA_T_BEGIN"] + [f"{g} {SIGMA_T[g]:.17g}" for g in range(2)] + ["SIGMA_T_END"]
    lines.append("TRANSFER_MOMENTS_BEGIN")
    for g in range(2):
        for gp in range(2):
            if S0[g][gp] != 0.0:
                lines.append(f"M_GFROM_GTO_VAL 0 {gp} {g} {S0[g][gp]:.17g}")
    lines.append("TRANSFER_MOMENTS_END")
    lines += ["SIGMA_F_BEGIN"] + [f"{g} {SIGMA_F[g]:.17g}" for g in range(2)] + ["SIGMA_F_END"]
    if scenario in ("matrix", "zero_column"):
        F = fission_operator(scenario)
        lines += ["CHI_BEGIN", "0 1.0", "1 0.0", "CHI_END", "PRODUCTION_MATRIX_BEGIN"]
        # Destination group first (see the user guide).
        lines += [f"GPRIME_G_VAL {g} {gp} {F[g][gp]:.17g}" for g in range(2) for gp in range(2)]
        lines.append("PRODUCTION_MATRIX_END")
    else:
        nu_prompt = NU_PROMPT_SIGMA_F / SIGMA_F
        nu_delayed = NU_DELAYED_SIGMA_F / SIGMA_F
        lines.insert(2, "NUM_PRECURSORS 1")
        lines += ["NU_PROMPT_BEGIN"] + [f"{g} {nu_prompt[g]:.17g}" for g in range(2)]
        lines += ["NU_PROMPT_END", "NU_DELAYED_BEGIN"]
        lines += [f"{g} {nu_delayed[g]:.17g}" for g in range(2)] + ["NU_DELAYED_END"]
        lines += ["CHI_PROMPT_BEGIN"] + [f"{g} {CHI_PROMPT[g]:.17g}" for g in range(2)]
        lines += ["CHI_PROMPT_END", "PRECURSOR_DECAY_CONSTANTS_BEGIN", "0 0.1"]
        lines += ["PRECURSOR_DECAY_CONSTANTS_END", "PRECURSOR_FRACTIONAL_YIELDS_BEGIN", "0 1.0"]
        lines += ["PRECURSOR_FRACTIONAL_YIELDS_END", "CHI_DELAYED_BEGIN"]
        lines += [f"G_PRECURSOR_VAL {g} 0 {CHI_DELAYED[g]:.17g}" for g in range(2)]
        lines += ["CHI_DELAYED_END", "INV_VELOCITY_BEGIN", "0 1.0", "1 1.0", "INV_VELOCITY_END"]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def remove_rank_file(prefix):
    try:
        os.remove(prefix + str(rank) + ".h5")
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    scenario = globals().get("scenario", "matrix")
    k_exact, production_exact, sigma_t_exact, scatter_exact = exact_relative_sensitivities(
        scenario
    )

    prefix = f"xs_sens_fission_op_{scenario}_"
    xs_file = prefix + "xs.xs"
    if rank == 0:
        write_xs(xs_file, scenario)
    barrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_file)
    barrier()
    if rank == 0:
        os.remove(xs_file)

    nodes = [i * 0.1 for i in range(11)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    boundaries = [{"name": "zmin", "type": "reflecting"}, {"name": "zmax", "type": "reflecting"}]
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 500,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=boundaries,
        options={
            "use_precursors": scenario == "delayed_on",
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
        **({"use_gpus": True} if globals().get("use_gpus", False) else {}),
    )

    def power_iteration():
        solver = PowerIterationKEigenSolver(problem=problem, max_iters=5000, k_tol=1.0e-12)
        solver.Initialize()
        solver.Execute()
        return solver.GetEigenvalue()

    fwd_phi, fwd_psi = prefix + "fwd_phi_", prefix + "fwd_psi_"
    adj_phi, adj_psi = prefix + "adj_phi_", prefix + "adj_psi_"
    k_eff = power_iteration()
    problem.WriteFluxMoments(fwd_phi)
    problem.WriteAngularFluxes(fwd_psi)
    problem.SetAdjoint(True)
    power_iteration()
    problem.WriteFluxMoments(adj_phi)
    problem.WriteAngularFluxes(adj_psi)

    def relative_k_sensitivity(**options):
        pp = CrossSectionSensitivityPostprocessor(
            problem=problem,
            relative=True,
            forward_flux_moments=fwd_phi,
            adjoint_flux_moments=adj_phi,
            forward_angular_fluxes=fwd_psi,
            adjoint_angular_fluxes=adj_psi,
            **options,
        )
        pp.Execute()
        pp.ApplyKEigenvalueScaling(k_eff)
        return pp.GetValue()[0][0] / k_eff

    def rel_err(value, exact):
        return abs(value - exact) / abs(exact)

    sigma_t_scale = max(abs(e) for e in sigma_t_exact)

    def sens_err(value, exact):
        # A sensitivity that is exactly 0, such as that of a group with no importance, has no
        # relative error, so it is measured against the largest sigma_t sensitivity.
        if abs(exact) <= 1.0e-12 * sigma_t_scale:
            return abs(value) / sigma_t_scale
        return rel_err(value, exact)

    production_err = max(
        rel_err(relative_k_sensitivity(sensitivity_type="production", group=gp),
                production_exact[gp])
        for gp in range(2)
        if production_exact[gp] != 0.0
    )
    sigma_t_err = max(
        sens_err(relative_k_sensitivity(sensitivity_type="sigma_t", group=g), sigma_t_exact[g])
        for g in range(2)
    )
    scatter_err = sens_err(
        relative_k_sensitivity(sensitivity_type="scatter", moment=0, from_group=0, to_group=1),
        scatter_exact,
    )

    def scaling_guarded():
        pp = CrossSectionSensitivityPostprocessor(
            problem=problem,
            sensitivity_type="production",
            group=0,
            forward_flux_moments=fwd_phi,
            adjoint_flux_moments=adj_phi,
        )

        def raises():
            try:
                pp.ApplyKEigenvalueScaling(k_eff)
            except RuntimeError:
                return True
            return False

        before_execute = raises()
        pp.Execute()
        pp.ApplyKEigenvalueScaling(k_eff)
        scaled = pp.GetValue()[0][0]
        twice = raises() and pp.GetValue()[0][0] == scaled
        pp.Execute()
        pp.ApplyKEigenvalueScaling(k_eff)
        return int(before_execute and twice and pp.GetValue()[0][0] == scaled)

    k_scaling_guarded = scaling_guarded()

    def zero_column_checks():
        relative = relative_k_sensitivity(sensitivity_type="production", group=1)
        pp = CrossSectionSensitivityPostprocessor(
            problem=problem,
            sensitivity_type="production",
            group=1,
            forward_flux_moments=fwd_phi,
            adjoint_flux_moments=adj_phi,
        )
        try:
            pp.Execute()
            rejected = False
        except ValueError as error:
            rejected = "undefined" in str(error)
        return relative, int(rejected)

    if scenario == "zero_column":
        zero_relative, zero_absolute_rejected = zero_column_checks()

    if rank == 0:
        print(f"XS_SENS_K_SCALING_GUARDED={k_scaling_guarded}")
        print(f"XS_SENS_FISSION_OP_K_REL_ERR={rel_err(k_eff, k_exact):.6e}")
        print(f"XS_SENS_FISSION_OP_PRODUCTION_REL_ERR={production_err:.6e}")
        print(f"XS_SENS_FISSION_OP_SIGMA_T_REL_ERR={sigma_t_err:.6e}")
        print(f"XS_SENS_FISSION_OP_SCATTER_REL_ERR={scatter_err:.6e}")
        if scenario == "zero_column":
            print(f"XS_SENS_ZERO_COLUMN_RELATIVE={zero_relative:.6e}")
            print(f"XS_SENS_ZERO_COLUMN_ABSOLUTE_REJECTED={zero_absolute_rejected}")

    barrier()
    for p in (fwd_phi, fwd_psi, adj_phi, adj_psi):
        remove_rank_file(p)
