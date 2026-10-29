#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Delayed neutrons in steady-state solves and in the transient that follows them.

Two groups with upscatter and two precursor families whose delayed spectra differ from each other
and from the prompt spectrum: chi_p = (0.9, 0.1), chi_1 = (0.5, 0.5), chi_2 = (0.1, 0.9), with
yields 0.4 and 0.6. A 1D slab with reflecting faces is an infinite medium: the flux is flat for any
quadrature and mesh, so numpy gives exact references. With L = diag(sigma_t) - S (S the transfer
matrix, row: to-group), F_p = chi_p (nu_p sigma_f)^T and F_d = sum_j gamma_j chi_j (nu_d sigma_f)^T,
use_precursors=True uses F = F_p + F_d and use_precursors=False uses F = F_p in every solve:

  steady_source: the fixed-source flux (L - F)^-1 q of a subcritical medium, with one and two
    groupsets.
  keigen: power iteration gives the dominant eigenvalue of L^-1 F: k_inf = 1.1262260 with delayed
    neutrons, k_p = 1.1051282 without. Its precursors, C_j = gamma_j nu_d sigma_f phi /
    (k lambda_j), are at equilibrium for the cross sections with nu_p and nu_d divided by k, which
    make the eigenstate critical. After that swap a transient must stay stationary for
    theta = 1 and 1/2 with either model.
  prompt_step: without the swap the prompt eigenstate is supercritical, and the first backward
    Euler step of the prompt-only transient solves (V^-1/dt + L - F_p) phi_1 = V^-1/dt phi_0.

Checks:
  DELAYED_INF_STEADY_SOURCE_MAX_REL_ERR: max over the four steady_source cases of
    max |phi - phi_ref| / max phi_ref.
  DELAYED_INF_KEIGEN_K_REL_ERR: max over the cases of |k / k_ref - 1|.
  DELAYED_INF_NULL_TRANSIENT_MAX_REL_DRIFT: max over models, theta, steps and nodes of
    |phi^n - phi^0| / max phi^0.
  DELAYED_INF_PROMPT_STEP_REL_ERR: max |phi_1 - phi_1,ref| / max phi_1,ref.
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return comm.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.source import VolumetricSource
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        SteadyStateSourceSolver,
        TransientSolver,
    )

NUM_GROUPS = 2
# (from-group, to-group, value), with upscatter
TRANSFERS = [(0, 0, 0.6), (0, 1, 0.25), (1, 1, 1.5), (1, 0, 0.02)]
MEDIUM = {
    "sigma_t": [1.0, 2.0],
    "sigma_f": [0.05, 0.25],
    "nu_prompt": [2.5, 2.4],
    "nu_delayed": [0.05, 0.04],
    "chi_prompt": [0.9, 0.1],
    "inv_velocity": [0.1, 1.0],
    # (emission spectrum, fractional yield, decay constant)
    "families": [([0.5, 0.5], 0.4, 0.0124), ([0.1, 0.9], 0.6, 0.305)],
}
DT = 0.1


def scaled(medium, sigma_f=1.0, nu=1.0):
    """The medium with sigma_f and with nu_p and nu_d multiplied by the given factors."""
    result = dict(medium)
    result["sigma_f"] = [sigma_f * v for v in medium["sigma_f"]]
    result["nu_prompt"] = [nu * v for v in medium["nu_prompt"]]
    result["nu_delayed"] = [nu * v for v in medium["nu_delayed"]]
    return result


def operators(medium):
    """L = diag(sigma_t) - S, F_p, F_d and V^-1."""
    S = np.zeros((NUM_GROUPS, NUM_GROUPS))
    for gf, gt, v in TRANSFERS:
        S[gt, gf] = v
    sigma_f = np.array(medium["sigma_f"])
    F_p = np.outer(medium["chi_prompt"], np.array(medium["nu_prompt"]) * sigma_f)
    nusf_d = np.array(medium["nu_delayed"]) * sigma_f
    F_d = sum(yld * np.outer(spectrum, nusf_d) for spectrum, yld, _ in medium["families"])
    return np.diag(medium["sigma_t"]) - S, F_p, F_d, np.diag(medium["inv_velocity"])


def k_infinite(medium, delayed=True):
    L, F_p, F_d, _ = operators(medium)
    values = np.linalg.eigvals(np.linalg.solve(L, F_p + F_d if delayed else F_p))
    return float(np.max(values.real))


def load_xs(medium, tag):
    """Writes the medium as an OpenSn cross-section file and loads it."""
    families = medium["families"]
    lines = [f"NUM_GROUPS {NUM_GROUPS}", "NUM_MOMENTS 1", f"NUM_PRECURSORS {len(families)}"]

    def block(name, values):
        lines.append(f"{name}_BEGIN")
        lines.extend(f"{g} {v!r}" for g, v in enumerate(values))
        lines.append(f"{name}_END")

    absorption = list(medium["sigma_t"])
    for gf, _, v in TRANSFERS:
        absorption[gf] -= v
    block("SIGMA_T", medium["sigma_t"])
    block("SIGMA_A", absorption)
    block("SIGMA_F", medium["sigma_f"])
    block("NU_PROMPT", medium["nu_prompt"])
    block("NU_DELAYED", medium["nu_delayed"])
    block("CHI_PROMPT", medium["chi_prompt"])
    block("INV_VELOCITY", medium["inv_velocity"])
    lines.append("CHI_DELAYED_BEGIN")
    for j, (spectrum, _, _) in enumerate(families):
        lines.extend(f"G_PRECURSOR_VAL {g} {j} {v!r}" for g, v in enumerate(spectrum))
    lines.append("CHI_DELAYED_END")
    block("PRECURSOR_DECAY_CONSTANTS", [decay for _, _, decay in families])
    block("PRECURSOR_FRACTIONAL_YIELDS", [yld for _, yld, _ in families])
    lines.append("TRANSFER_MOMENTS_BEGIN")
    lines.extend(f"M_GFROM_GTO_VAL 0 {gf} {gt} {v!r}" for gf, gt, v in TRANSFERS)
    lines.append("TRANSFER_MOMENTS_END")

    path = f"transient_delayed_production_infinite_medium_{tag}_r{rank}.xs"
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    os.remove(path)
    return xs


def build_problem(xs, use_precursors, groupsets=1, source=None):
    nodes = [0.25 * i for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    group_ranges = [[0, 1]] if groupsets == 1 else [[0, 0], [1, 1]]
    options = {
        "save_angular_flux": True,
        "use_precursors": use_precursors,
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": False,
    }
    if groupsets > 1:
        options["ags_tolerance"] = 1.0e-13
        options["max_ags_iterations"] = 1000
    settings = {}
    if source is not None:
        settings["volumetric_sources"] = [VolumetricSource(block_ids=[0], group_strength=source)]
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=NUM_GROUPS,
        groupsets=[
            {
                "groups_from_to": group_range,
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 200,
            }
            for group_range in group_ranges
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options=options,
        **settings,
    )


def group_flux(problem):
    """The scalar flux by node and group."""
    return np.array(problem.GetPhiNewLocal()).reshape(-1, NUM_GROUPS)


def max_rel_err(phi, phi_ref):
    """max |phi - phi_ref| / max phi_ref over all ranks, nodes and groups."""
    local = float(np.max(np.abs(phi - phi_ref))) if phi.size else 0.0
    return MPIAllReduce(local, "max") / float(np.max(np.abs(phi_ref)))


def keigen_state(use_precursors, tag):
    """A problem with its PI eigenstate made critical by dividing nu by k, in transient mode."""
    problem = build_problem(load_xs(MEDIUM, tag), use_precursors)
    pi = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-13, max_iters=1000)
    pi.Initialize()
    pi.Execute()
    k = pi.GetEigenvalue()
    critical_xs = load_xs(scaled(MEDIUM, nu=1.0 / k), tag + "_critical")
    return problem, k, critical_xs


if __name__ == "__main__":
    results = {}

    # Fixed source, subcritical (k_inf = 0.68)
    source_medium = scaled(MEDIUM, sigma_f=0.6)
    q = np.array([1.0, 0.5])
    L, F_p, F_d, _ = operators(source_medium)
    err = 0.0
    for use_precursors in (True, False):
        phi_ref = np.linalg.solve(L - F_p - (F_d if use_precursors else 0.0), q)
        for groupsets in (1, 2):
            problem = build_problem(
                load_xs(source_medium, f"source_{use_precursors}_{groupsets}"),
                use_precursors,
                groupsets,
                source=q.tolist(),
            )
            solver = SteadyStateSourceSolver(problem=problem)
            solver.Initialize()
            solver.Execute()
            phi = group_flux(problem)
            err = max(err, max_rel_err(phi, np.broadcast_to(phi_ref, phi.shape)))
    results["STEADY_SOURCE_MAX_REL_ERR"] = err

    # k-eigenvalue and null transients from the critical eigenstate
    k_err = 0.0
    drift = 0.0
    for use_precursors in (True, False):
        k_ref = k_infinite(MEDIUM, delayed=use_precursors)
        for theta in (1.0, 0.5):
            tag = f"keigen_{use_precursors}_{theta}"
            problem, k, critical_xs = keigen_state(use_precursors, tag)
            k_err = max(k_err, abs(k / k_ref - 1.0))
            phi0 = group_flux(problem).copy()
            problem.SetTimeDependentMode()
            transient = TransientSolver(
                problem=problem, dt=DT, stop_time=1.0e9, theta=theta, initial_state="existing",
                verbose=False,
            )
            transient.Initialize()
            problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": critical_xs}])
            for _ in range(5):
                transient.Advance()
                drift = max(drift, max_rel_err(group_flux(problem), phi0))
    results["KEIGEN_K_REL_ERR"] = k_err
    results["NULL_TRANSIENT_MAX_REL_DRIFT"] = drift

    # Prompt-only transient step from the supercritical prompt eigenstate
    problem, k, _ = keigen_state(False, "prompt_step")
    phi0 = group_flux(problem).copy()
    problem.SetTimeDependentMode()
    transient = TransientSolver(
        problem=problem, dt=DT, stop_time=1.0e9, theta=1.0, initial_state="existing",
        verbose=False,
    )
    transient.Initialize()
    transient.Advance()
    phi1 = group_flux(problem)
    L, F_p, _, V_inv = operators(MEDIUM)
    # phi_0 is flat: every node has the same group flux.
    phi0_group = np.array(
        [MPIAllReduce(float(v), "max") for v in (phi0[0] if phi0.size else np.zeros(NUM_GROUPS))]
    )
    phi1_ref = np.linalg.solve(V_inv / DT + L - F_p, V_inv @ phi0_group / DT)
    results["PROMPT_STEP_REL_ERR"] = max_rel_err(phi1, np.broadcast_to(phi1_ref, phi1.shape))

    if rank == 0:
        for key, value in results.items():
            print(f"DELAYED_INF_{key}={value:.10e}")
