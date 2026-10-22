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
  no_precursors: use_precursors=False, so no delayed neutrons are produced: k_inf and the P/R
    references are those of the prompt-only F = chi_p (nu_p sigma_f)^T (k_inf = 0.993), and k and
    P/R are checked as in production.

Two scenarios use two materials, A and B, which alternate cell by cell on a mesh of 20 cells. They
have the same sigma_t, S, sigma_f and nu = nu_p + nu_d, and the same F, but split F differently
into its prompt and delayed parts: A has delayed fractions 0.006 and 0.008 (of nu in the fissioning
group) and chi_d = (0.675, 0.325), B about 0.025 and chi_d = (0.1, 0.9), each from two families of
its own. B's prompt spectrum is (0.98, 0.02) against A's (0.96, 0.04), and B's delayed fractions
are those for which F_B = F_A. The slab is then an infinite medium of F, with k_inf = 0.8701171 and
a flat flux, provided that each cell uses the delayed spectrum of its own material. With A's
spectrum in B's cells, F there changes by 1.5e-2 of its largest entry (by 4.8e-3 in A's cells with
B's spectrum).
  two_materials: forward and adjoint k and P/R as in production, and the flux must be flat: the
    largest (max - min) / max |phi| over the groups, the extremes taken over the nodes and
    quadrature points of the slab.
  two_materials_source: the same slab as a fixed-source problem with one groupset per group,
    a flat source q = (1, 0.25), and the flat adjoint source q* = (0.3, 1). With the precursors
    in equilibrium the fluxes are flat: phi = (L - F)^-1 q and phi* = (L - F)^-T q*, with
    L = diag(sigma_t) - S. Unlike power iteration, this solve evaluates the fission source
    within a groupset and across groupsets in separate calls. Checked are the fluxes against
    phi and phi* and the ratio P/R of each against 1^T F phi / sigma_f^T phi and
    1^T F^T phi* / sigma_f^T phi*.

Checks:
  ADJOINT_DELAYED_K_FWD_REL_ERR: |k_forward / k_inf - 1|.
  ADJOINT_DELAYED_K_ADJ_REL_ERR: |k_adjoint / k_inf - 1|.
  ADJOINT_DELAYED_<SCENARIO>_K_{FWD,ADJ}_REL_ERR: the same relative errors for the other
    k-eigenvalue scenarios.
  ADJOINT_DELAYED_<SCENARIO>_RATIO_{FWD,ADJ}_REL_ERR: the relative errors of P/R (production,
    no_precursors, two_materials and two_materials_source).
  ADJOINT_DELAYED_TWO_MATERIALS_FLAT_{FWD,ADJ}: the flatness of the forward and adjoint flux.
  ADJOINT_DELAYED_TWO_MATERIALS_SOURCE_PHI_{FWD,ADJ}_REL_ERR: the largest |phi - phi_ref| over
    the slab and the groups, divided by the largest phi_ref.
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
    from pyopensn.solver import CMFDAcceleration, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

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


def split_medium(chi_prompt, families, beta, nu=(2.0, 2.2)):
    """A medium of the two-material scenarios: delayed fraction beta[g'] of nu[g']."""
    return {
        "sigma_t": [1.0, 2.0],
        "sigma_a": [0.2, 0.5],
        "sigma_f": [0.05, 0.25],
        "nu_prompt": [(1.0 - b) * n for b, n in zip(beta, nu)],
        "nu_delayed": [b * n for b, n in zip(beta, nu)],
        "chi_prompt": chi_prompt,
        "families": families,
    }


def delayed_spectrum(medium):
    """chi_d = sum_j gamma_j chi_j."""
    return sum(yld * np.array(spectrum) for spectrum, yld, _ in medium["families"])


def two_materials():
    """Materials A and B of the two-material scenarios, with the same F."""
    a_beta = np.array([0.006, 0.008])
    a = split_medium(
        [0.96, 0.04], [([0.9, 0.1], 0.25, 0.0127), ([0.6, 0.4], 0.75, 0.317)], a_beta.tolist()
    )
    b_chi_prompt = [0.98, 0.02]
    b_families = [([0.05, 0.95], 0.5, 0.0335), ([0.15, 0.85], 0.5, 1.22)]
    # A neutron from fission in group g' is born fast with the same probability in both:
    # chi_p[0] + beta[g'] (chi_d[0] - chi_p[0]). This gives B's delayed fractions.
    a_fast = a["chi_prompt"][0] + a_beta * (delayed_spectrum(a)[0] - a["chi_prompt"][0])
    b_chi_d = delayed_spectrum({"families": b_families})
    b_beta = (b_chi_prompt[0] - a_fast) / (b_chi_prompt[0] - b_chi_d[0])
    b = split_medium(b_chi_prompt, b_families, b_beta.tolist())
    return [a, b]


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


def production_ratio(medium, phi, adjoint, use_precursors=True):
    """1^T F phi / sigma_f^T phi, or with F^T in adjoint mode."""
    _, F = infinite_medium_operators(medium, use_precursors)
    if adjoint:
        F = F.T
    return float(np.sum(F @ phi) / (np.array(medium["sigma_f"]) @ phi))


def production_ratio_reference(medium, adjoint, use_precursors=True):
    """The production ratio of the forward eigenvector, or of the adjoint one."""
    L, F = infinite_medium_operators(medium, use_precursors)
    if adjoint:
        L, F = L.T, F.T
    _, phi = dominant(L, F)
    return production_ratio(medium, phi, adjoint, use_precursors)


def flux_range(phys):
    """The smallest and largest scalar flux of each group, over the nodes and quadrature points."""
    everywhere = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ranges = []
    for ff in phys.GetScalarFluxFieldFunction(only_scalar_flux=True):
        values = []
        for operation in ("min", "max"):
            ffi = FieldFunctionInterpolationVolume()
            ffi.SetOperationType(operation)
            ffi.SetLogicalVolume(everywhere)
            ffi.AddFieldFunction(ff)
            ffi.Execute()
            values.append(ffi.GetValue())
        ranges.append(values)
    return np.array(ranges)


def solve(grid, materials, adjoint, groupsets=1, cmfd=False, use_precursors=True, source=None):
    """Material i on block i. With a source (group strengths), a fixed-source solve."""
    # Each problem has its own cross-section objects, and each scenario its own files.
    xs_map = []
    for block_id, medium in enumerate(materials):
        path = f"adjoint_keigen_delayed_infinite_medium_{scenario}_{block_id}_r{rank}.xs"
        write_xs(path, medium)
        xs = MultiGroupXS()
        xs.LoadFromOpenSn(path)
        os.remove(path)
        xs_map.append({"block_ids": [block_id], "xs": xs})

    group_ranges = [[0, 1]] if groupsets == 1 else [[0, 0], [1, 1]]
    options = {
        "adjoint": adjoint,
        "use_precursors": use_precursors,
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": False,
    }
    if groupsets > 1:
        options["ags_tolerance"] = 1.0e-12
    problem_settings = {}
    if source is not None:
        options["max_ags_iterations"] = 1000
        problem_settings["volumetric_sources"] = [
            VolumetricSource(block_ids=list(range(len(materials))), group_strength=source)
        ]
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
        xs_map=xs_map,
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options=options,
        **problem_settings,
    )
    if source is not None:
        solver = SteadyStateSourceSolver(problem=phys)
    elif cmfd:
        # The settings of keigenvalue_transport_1d_1g_cmfd.py.
        acceleration = CMFDAcceleration(
            problem=phys,
            current_closure="partial",
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
    k = None if source is not None else solver.GetEigenvalue()
    ratio = phys.ComputeFissionProduction("new") / phys.ComputeFissionRate("new")
    return k, ratio, flux_range(phys)


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
        k_forward, _, _ = solve(grid, [medium], adjoint=False)
        k_adjoint, _, _ = solve(grid, [medium], adjoint=True)
        prefix = "ADJOINT_DELAYED_"
        results["K_INF="] = f"{k_inf:.10f}"
        results["K_FWD="] = f"{k_forward:.10f}"
        results["K_ADJ="] = f"{k_adjoint:.10f}"
        results["K_FWD_REL_ERR="] = f"{abs(k_forward / k_inf - 1.0):.6e}"
        results["K_ADJ_REL_ERR="] = f"{abs(k_adjoint / k_inf - 1.0):.6e}"
    elif scenario in ("two_groupsets", "cmfd", "production", "no_precursors"):
        medium = media["two_families"]
        # Without precursors the references are those of the prompt-only operator.
        use_precursors = scenario != "no_precursors"
        k_inf = k_infinite(medium, use_precursors)
        settings = {
            "groupsets": 2 if scenario == "two_groupsets" else 1,
            "cmfd": scenario == "cmfd",
            "use_precursors": use_precursors,
        }
        k_forward, ratio_forward, _ = solve(grid, [medium], adjoint=False, **settings)
        k_adjoint, ratio_adjoint, _ = solve(grid, [medium], adjoint=True, **settings)
        prefix = f"ADJOINT_DELAYED_{scenario.upper()}_"
        results["K_INF="] = f"{k_inf:.10f}"
        results["K_FWD="] = f"{k_forward:.10f}"
        results["K_ADJ="] = f"{k_adjoint:.10f}"
        results["K_FWD_REL_ERR="] = f"{abs(k_forward / k_inf - 1.0):.6e}"
        results["K_ADJ_REL_ERR="] = f"{abs(k_adjoint / k_inf - 1.0):.6e}"
        if scenario in ("production", "no_precursors"):
            ratio_forward_ref = production_ratio_reference(
                medium, adjoint=False, use_precursors=use_precursors
            )
            ratio_adjoint_ref = production_ratio_reference(
                medium, adjoint=True, use_precursors=use_precursors
            )
            results["RATIO_FWD="] = f"{ratio_forward:.10f}"
            results["RATIO_ADJ="] = f"{ratio_adjoint:.10f}"
            results["RATIO_FWD_REF="] = f"{ratio_forward_ref:.10f}"
            results["RATIO_ADJ_REF="] = f"{ratio_adjoint_ref:.10f}"
            results["RATIO_FWD_REL_ERR="] = f"{abs(ratio_forward / ratio_forward_ref - 1.0):.6e}"
            results["RATIO_ADJ_REL_ERR="] = f"{abs(ratio_adjoint / ratio_adjoint_ref - 1.0):.6e}"
    elif scenario in ("two_materials", "two_materials_source"):
        materials = two_materials()
        L, F = infinite_medium_operators(materials[0])
        F_b = infinite_medium_operators(materials[1])[1]
        if np.abs(F_b - F).max() > 1.0e-14 * np.abs(F).max():
            sys.exit("The two materials must have the same production matrix.")
        medium = materials[0]

        # 20 cells, 0.05 wide, alternately in A (block 0) and B (block 1)
        nodes = [0.05 * i for i in range(21)]
        grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
        grid.SetBlockIDFromFunction(lambda centroid, block_id: int(centroid.z / 0.05) % 2)
        prefix = f"ADJOINT_DELAYED_{scenario.upper()}_"
        if scenario == "two_materials":
            k_inf = k_infinite(medium)
            k_forward, ratio_forward, range_forward = solve(grid, materials, adjoint=False)
            k_adjoint, ratio_adjoint, range_adjoint = solve(grid, materials, adjoint=True)
            ratio_forward_ref = production_ratio_reference(medium, adjoint=False)
            ratio_adjoint_ref = production_ratio_reference(medium, adjoint=True)
            results["K_INF="] = f"{k_inf:.10f}"
            results["K_FWD="] = f"{k_forward:.10f}"
            results["K_ADJ="] = f"{k_adjoint:.10f}"
            results["K_FWD_REL_ERR="] = f"{abs(k_forward / k_inf - 1.0):.6e}"
            results["K_ADJ_REL_ERR="] = f"{abs(k_adjoint / k_inf - 1.0):.6e}"
            for name, flux in (("FWD", range_forward), ("ADJ", range_adjoint)):
                flatness = np.max((flux[:, 1] - flux[:, 0]) / np.abs(flux).max(axis=1))
                results[f"FLAT_{name}="] = f"{flatness:.6e}"
        else:
            q_forward = [1.0, 0.25]
            q_adjoint = [0.3, 1.0]
            _, ratio_forward, range_forward = solve(
                grid, materials, adjoint=False, groupsets=2, source=q_forward
            )
            _, ratio_adjoint, range_adjoint = solve(
                grid, materials, adjoint=True, groupsets=2, source=q_adjoint
            )
            phi_forward = np.linalg.solve(L - F, q_forward)
            phi_adjoint = np.linalg.solve((L - F).T, q_adjoint)
            ratio_forward_ref = production_ratio(medium, phi_forward, adjoint=False)
            ratio_adjoint_ref = production_ratio(medium, phi_adjoint, adjoint=True)
            for name, flux, ref in (
                ("FWD", range_forward, phi_forward),
                ("ADJ", range_adjoint, phi_adjoint),
            ):
                error = np.abs(flux - ref[:, np.newaxis]).max() / np.abs(ref).max()
                results[f"PHI_{name}_MIN="] = " ".join(f"{v:.10f}" for v in flux[:, 0])
                results[f"PHI_{name}_MAX="] = " ".join(f"{v:.10f}" for v in flux[:, 1])
                results[f"PHI_{name}_REF="] = " ".join(f"{v:.10f}" for v in ref)
                results[f"PHI_{name}_REL_ERR="] = f"{error:.6e}"
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
