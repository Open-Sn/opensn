#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Adjoint k-eigenvalue with delayed-neutron precursors in an infinite medium.

Two groups, one precursor family, prompt neutrons born fast (chi_p = (1, 0)) and delayed neutrons
born in both groups (chi_d = (0.2, 0.8)), delayed fraction 0.0065, use_precursors=True. A 1D slab
with both faces reflecting is an infinite medium: the flat flux is exact for any quadrature and
mesh, so k is the dominant eigenvalue k_inf of (diag(sigma_t) - S)^-1 F, with S the transfer matrix
(row: to-group) and F = chi_p (nu_p sigma_f)^T + chi_d (nu_d sigma_f)^T. The adjoint problem
transposes S and F, and its eigenvalue is the same k_inf. With the delayed part of F left
untransposed, the adjoint eigenvalue would be 0.9676877 instead of k_inf = 0.967488.

The other scenarios (set `scenario`) use a second medium: two precursor families with different
emission spectra and yields 0.4 and 0.6, group-dependent nu, and prompt neutrons born in both
groups, so that the family, its yield and both group indices of the delayed term matter
(k_inf = 1.0005758; untransposed delayed part: 1.0007892).
  two_groupsets: groupsets [0] and [1]; the delayed term between the groups enters through the
    across-groupset source.
  cmfd: power iteration with CMFD acceleration, whose coarse production includes the delayed
    term. On this slab CMFD converges to about 1e-6 of k_inf, so the tolerance is 1e-5; an
    untransposed coarse delayed term moves the eigenvalue by 2.1e-4.
  production: the ratio P/R of the fission production (ComputeFissionProduction) to the fission
    rate (ComputeFissionRate) of the converged flux, against numpy's 1^T F phi / sigma_f^T phi
    (forward) and 1^T F^T phi* / sigma_f^T phi* (adjoint), phi* the eigenvector of the transposed
    operator. k does not depend on how the production weights the flux; the ratio does.
  no_precursors: use_precursors=False; the adjoint eigenvalue equals the forward one. The forward
    eigenvalue selects the operator, F with or without the delayed part, whichever k is nearer, so
    the check does not fix whether such a solve applies delayed production. P/R of the forward and
    adjoint fluxes is compared with that operator's references, as in production.

Checks:
  ADJOINT_DELAYED_K_FWD_REL_ERR: |k_forward / k_inf - 1|.
  ADJOINT_DELAYED_K_ADJ_REL_ERR: |k_adjoint / k_inf - 1|.
  ADJOINT_DELAYED_<SCENARIO>_K_{FWD,ADJ}_REL_ERR: the same relative errors for the scenarios
    two_groupsets, cmfd and production.
  ADJOINT_DELAYED_PRODUCTION_RATIO_{FWD,ADJ}_REL_ERR: the relative errors of P/R.
  ADJOINT_DELAYED_NO_PRECURSORS_K_ADJ_FWD_REL_DIFF: |k_adjoint / k_forward - 1|.
  ADJOINT_DELAYED_NO_PRECURSORS_K_FWD_REL_ERR: |k_forward / k_ref - 1|, k_ref the selected
    operator's eigenvalue.
  ADJOINT_DELAYED_NO_PRECURSORS_RATIO_{FWD,ADJ}_REL_ERR: the relative errors of P/R against the
    selected operator.
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration

# (from-group, to-group, value)
transfers = [(0, 0, 0.5), (0, 1, 0.3), (1, 1, 1.5)]
media = {
    "one_family": {
        "sigma_t": [1.0, 2.0],
        "sigma_a": [0.2, 0.5],
        "sigma_f": [0.05, 0.25],
        "nu_prompt": [2.4, 2.4],
        "nu_delayed": [0.0156, 0.0156],
        "chi_prompt": [1.0, 0.0],
        # (emission spectrum, fractional yield, decay constant)
        "families": [([0.2, 0.8], 1.0, 0.08)],
    },
    "two_families": {
        "sigma_t": [1.0, 2.0],
        "sigma_a": [0.2, 0.5],
        "sigma_f": [0.05, 0.25],
        "nu_prompt": [2.5, 2.4],
        "nu_delayed": [0.02, 0.015],
        "chi_prompt": [0.9, 0.1],
        "families": [([0.5, 0.5], 0.4, 0.0124), ([0.1, 0.9], 0.6, 0.305)],
    },
}


def write_xs(path, medium):
    families = medium["families"]
    lines = ["NUM_GROUPS 2", "NUM_MOMENTS 1", f"NUM_PRECURSORS {len(families)}"]

    def block(name, values):
        lines.append(f"{name}_BEGIN")
        lines.extend(f"{g} {v!r}" for g, v in enumerate(values))
        lines.append(f"{name}_END")

    block("SIGMA_T", medium["sigma_t"])
    block("SIGMA_A", medium["sigma_a"])
    block("SIGMA_F", medium["sigma_f"])
    block("NU_PROMPT", medium["nu_prompt"])
    block("NU_DELAYED", medium["nu_delayed"])
    block("CHI_PROMPT", medium["chi_prompt"])
    lines.append("CHI_DELAYED_BEGIN")
    for j, (spectrum, _, _) in enumerate(families):
        lines.extend(f"G_PRECURSOR_VAL {g} {j} {v!r}" for g, v in enumerate(spectrum))
    lines.append("CHI_DELAYED_END")
    block("PRECURSOR_DECAY_CONSTANTS", [decay for _, _, decay in families])
    block("PRECURSOR_FRACTIONAL_YIELDS", [yld for _, yld, _ in families])
    lines.append("TRANSFER_MOMENTS_BEGIN")
    lines.extend(f"M_GFROM_GTO_VAL 0 {gf} {gt} {v!r}" for gf, gt, v in transfers)
    lines.append("TRANSFER_MOMENTS_END")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def infinite_medium_operators(medium, use_precursors=True):
    """Return L = diag(sigma_t) - S and the production matrix F (row: to-group)."""
    S = np.zeros((2, 2))
    for gf, gt, v in transfers:
        S[gt, gf] = v
    sigma_f = np.array(medium["sigma_f"])
    F = np.outer(medium["chi_prompt"], np.array(medium["nu_prompt"]) * sigma_f)
    if use_precursors:
        nusf_d = np.array(medium["nu_delayed"]) * sigma_f
        for spectrum, yld, _ in medium["families"]:
            F += yld * np.outer(spectrum, nusf_d)
    return np.diag(medium["sigma_t"]) - S, F


def dominant(L, F):
    """Dominant eigenvalue and eigenvector of L^-1 F."""
    values, vectors = np.linalg.eig(np.linalg.solve(L, F))
    i = np.argmax(values.real)
    return float(values[i].real), vectors[:, i].real


def k_infinite(medium, use_precursors=True):
    return dominant(*infinite_medium_operators(medium, use_precursors))[0]


def production_ratio_reference(medium, adjoint, use_precursors=True):
    """1^T F phi / sigma_f^T phi for the forward eigenvector, or with F^T for the adjoint one."""
    L, F = infinite_medium_operators(medium, use_precursors)
    if adjoint:
        L, F = L.T, F.T
    _, phi = dominant(L, F)
    return float(np.sum(F @ phi) / (np.array(medium["sigma_f"]) @ phi))


def solve(grid, medium, adjoint, groupsets=1, cmfd=False, use_precursors=True):
    # Each problem has its own cross-section object, and each scenario its own file.
    path = f"adjoint_keigen_delayed_infinite_medium_{scenario}_r{rank}.xs"
    write_xs(path, medium)
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    os.remove(path)

    group_ranges = [[0, 1]] if groupsets == 1 else [[0, 0], [1, 1]]
    options = {
        "adjoint": adjoint,
        "use_precursors": use_precursors,
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": False,
    }
    if groupsets > 1:
        options["ags_tolerance"] = 1.0e-12
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": group_range,
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
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
    )
    if cmfd:
        # The settings of keigenvalue_transport_1d_1g_cmfd.py.
        acceleration = CMFDAcceleration(
            problem=phys,
            current_closure="auto",
            coarse_mesh="local_aggregation",
            aggregation_size=2,
            group_aggregation_size=1,
            relaxation=0.5,
            inactive_iterations=0,
            update_wgs_max_its=1,
            update_wgs_abs_tol=1.0e-12,
            correction_max_attempts=10,
            correction_min_damping=1.0e-4,
            negative_flux_tolerance=1.0e-6,
            petsc_options="-CMFDAccelerationksp_type gmres -CMFDAccelerationpc_type jacobi",
        )
        solver = PowerIterationKEigenSolver(
            problem=phys, acceleration=acceleration, max_iters=2000, k_tol=1.0e-10
        )
    else:
        solver = PowerIterationKEigenSolver(problem=phys, k_tol=1.0e-12)
    solver.Initialize()
    solver.Execute()
    ratio = phys.ComputeFissionProduction("new") / phys.ComputeFissionRate("new")
    return solver.GetEigenvalue(), ratio


if __name__ == "__main__":
    num_procs = 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    if "scenario" not in globals():
        scenario = "one_groupset"

    nodes = [0.1 * i for i in range(11)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    results = {}
    if scenario == "one_groupset":
        medium = media["one_family"]
        k_inf = k_infinite(medium)
        k_forward, _ = solve(grid, medium, adjoint=False)
        k_adjoint, _ = solve(grid, medium, adjoint=True)
        prefix = "ADJOINT_DELAYED_"
        results["K_INF="] = f"{k_inf:.10f}"
        results["K_FWD="] = f"{k_forward:.10f}"
        results["K_ADJ="] = f"{k_adjoint:.10f}"
        results["K_FWD_REL_ERR="] = f"{abs(k_forward / k_inf - 1.0):.6e}"
        results["K_ADJ_REL_ERR="] = f"{abs(k_adjoint / k_inf - 1.0):.6e}"
    elif scenario in ("two_groupsets", "cmfd", "production"):
        medium = media["two_families"]
        k_inf = k_infinite(medium)
        settings = {
            "groupsets": 2 if scenario == "two_groupsets" else 1,
            "cmfd": scenario == "cmfd",
        }
        k_forward, ratio_forward = solve(grid, medium, adjoint=False, **settings)
        k_adjoint, ratio_adjoint = solve(grid, medium, adjoint=True, **settings)
        prefix = f"ADJOINT_DELAYED_{scenario.upper()}_"
        results["K_INF="] = f"{k_inf:.10f}"
        results["K_FWD="] = f"{k_forward:.10f}"
        results["K_ADJ="] = f"{k_adjoint:.10f}"
        results["K_FWD_REL_ERR="] = f"{abs(k_forward / k_inf - 1.0):.6e}"
        results["K_ADJ_REL_ERR="] = f"{abs(k_adjoint / k_inf - 1.0):.6e}"
        if scenario == "production":
            ratio_forward_ref = production_ratio_reference(medium, adjoint=False)
            ratio_adjoint_ref = production_ratio_reference(medium, adjoint=True)
            results["RATIO_FWD="] = f"{ratio_forward:.10f}"
            results["RATIO_ADJ="] = f"{ratio_adjoint:.10f}"
            results["RATIO_FWD_REF="] = f"{ratio_forward_ref:.10f}"
            results["RATIO_ADJ_REF="] = f"{ratio_adjoint_ref:.10f}"
            results["RATIO_FWD_REL_ERR="] = f"{abs(ratio_forward / ratio_forward_ref - 1.0):.6e}"
            results["RATIO_ADJ_REL_ERR="] = f"{abs(ratio_adjoint / ratio_adjoint_ref - 1.0):.6e}"
    elif scenario == "no_precursors":
        medium = media["two_families"]
        k_forward, ratio_forward = solve(grid, medium, adjoint=False, use_precursors=False)
        k_adjoint, ratio_adjoint = solve(grid, medium, adjoint=True, use_precursors=False)
        k_prompt = k_infinite(medium, use_precursors=False)
        k_full = k_infinite(medium)
        delayed = abs(k_forward / k_full - 1.0) < abs(k_forward / k_prompt - 1.0)
        k_ref = k_full if delayed else k_prompt
        ratio_forward_ref = production_ratio_reference(
            medium, adjoint=False, use_precursors=delayed
        )
        ratio_adjoint_ref = production_ratio_reference(medium, adjoint=True, use_precursors=delayed)
        prefix = "ADJOINT_DELAYED_NO_PRECURSORS_"
        results["DELAYED_PRODUCTION="] = f"{delayed}"
        results["K_REF="] = f"{k_ref:.10f}"
        results["K_FWD="] = f"{k_forward:.10f}"
        results["K_ADJ="] = f"{k_adjoint:.10f}"
        results["K_ADJ_FWD_REL_DIFF="] = f"{abs(k_adjoint / k_forward - 1.0):.6e}"
        results["K_FWD_REL_ERR="] = f"{abs(k_forward / k_ref - 1.0):.6e}"
        results["RATIO_FWD="] = f"{ratio_forward:.10f}"
        results["RATIO_ADJ="] = f"{ratio_adjoint:.10f}"
        results["RATIO_FWD_REF="] = f"{ratio_forward_ref:.10f}"
        results["RATIO_ADJ_REF="] = f"{ratio_adjoint_ref:.10f}"
        results["RATIO_FWD_REL_ERR="] = f"{abs(ratio_forward / ratio_forward_ref - 1.0):.6e}"
        results["RATIO_ADJ_REL_ERR="] = f"{abs(ratio_adjoint / ratio_adjoint_ref - 1.0):.6e}"
    else:
        sys.exit(f"Unknown scenario {scenario!r}.")

    if rank == 0:
        for key, value in results.items():
            print(f"{prefix}{key}{value}")
