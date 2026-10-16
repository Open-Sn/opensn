#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Forward and adjoint solves on a cross-section object shared between problems or rescaled
between mode switches.

Two groups, sigma_t = (1, 1), within-group scattering 0.5 in each group, downscattering 0.3 from
group 0 to group 1. A 1D slab with both faces reflecting is an infinite medium: the flat flux is
exact for any quadrature and mesh, phi = (diag(sigma_t) - S)^-1 q in forward mode and
phi = (diag(sigma_t) - S^T)^-1 q in adjoint mode, with S the transfer matrix (row: to-group).

A. One cross-section object is used by a forward problem with q = (1, 0), phi = (2, 1.2), and then
   by an adjoint problem with q = (0, 1), phi = (1.2, 2). Solving the forward problem again must
   give (2, 1.2) (with transposed scattering it would give (2, 0)), and the adjoint problem
   (1.2, 2).
B. A problem is switched to adjoint mode and back, the cross sections are scaled by 2, and the
   problem is switched to adjoint mode again. Its adjoint solve with q = (0, 1) must use the scaled
   data, phi = (0.6, 1) (with the unscaled transfers it would give (0.1333, 0.6667)).
C. An adjoint problem with q = (0, 1) is built, then a forward problem with q = (1, 0) is built on
   the same cross sections, the cross sections are scaled by 2, and both problems are solved
   without being rebuilt. The adjoint solve must give (0.6, 1) (with the unscaled transposes,
   (0.1333, 0.6667)) and the forward solve (1, 0.6).

Checks: max over nodes and groups of |phi / phi_exact - 1| for
  SHARED_XS_FWD_FIRST_ERR, SHARED_XS_ADJ_ERR, SHARED_XS_FWD_AGAIN_ERR, SHARED_XS_ADJ_AGAIN_ERR (A),
  SCALED_XS_ADJ_ERR (B), and SHARED_XS_LIVE_ADJ_SCALED_ERR, SHARED_XS_LIVE_FWD_SCALED_ERR (C).

The other scenarios (set `scenario`) add fission to the same medium: sigma_f = (0.05, 0.25),
nu = (2.5, 2.4), chi = (0.9, 0.1), so F = chi (nu sigma_f)^T. The k-eigenvalue of both modes is
k_inf = 0.993, the dominant eigenvalue of (diag(sigma_t) - S)^-1 F; the forward operator with F^T
in place of F would give 0.36. The ratio P/R of the fission production (ComputeFissionProduction)
to the fission rate (ComputeFissionRate) of the converged flux must equal numpy's
1^T F phi / sigma_f^T phi (forward) or 1^T F^T phi* / sigma_f^T phi* (adjoint, phi* the
eigenvector of the transposed operator), 2.4220 and 2.0658; k does not depend on how the
production weights the flux, the ratio does.
  keigen: one cross-section object is used by a forward and an adjoint problem, built in either
    order (FA: forward first, AF: adjoint first); each problem is solved after both are built.
  keigen_cmfd: the FA order, with CMFD acceleration on the forward problem. On this slab the
    CMFD eigenvalue is about 1.5e-6 off k_inf, so the tolerance of that problem is 1e-5. Its power
    iteration must also stop before its iteration limit: when CMFD skips its corrections, the
    unaccelerated iteration reaches the limit with k_inf to roundoff.
  keigen_cmfd_groups: keigen_cmfd with group_aggregation_size = 2. Both groups form one coarse
    group, so the condensed removal and production cross sections contain the transfer between
    the two fine groups; with one fine group per coarse group only their diagonals enter.
  keigen_cmfd_adjoint: CMFD on the adjoint problem instead, with group_aggregation_size = 1 and 2,
    alone (A) and with a forward problem built later on the same cross sections (AF).
  keigen_delayed: the keigen orders with two delayed-neutron precursor families: nu = nu_p + nu_d
    with nu_p = (2.4, 2.3) and nu_d = (0.1, 0.1), chi_p = (0.9, 0.1), family spectra (0.3, 0.7) and
    (0.6, 0.4) with fractional yields (0.4, 0.6), so F = chi_p (nu_p sigma_f)^T + chi_d
    (nu_d sigma_f)^T with chi_d = (0.48, 0.52) and k_inf = 0.99699. The adjoint problem is built on
    the same cross sections, and only the forward problem is solved and checked; a transposed
    delayed term would put its k 7.7e-3 off.
  first_collision: on a 2D mesh with vacuum faces and the medium without fission, an adjoint
    problem is built, then a forward problem whose first-collision source comes from an
    uncollided-flux file, on the same cross sections, and both are solved. Each flux must equal
    that of the same problem on its own copy of the cross sections.
  first_collision_fission: first_collision with fission, so the first-collision source also
    contains the fission production of the uncollided flux.
  reload: the fixed sources of A, on the medium without fission and on the medium with sigma_f
    halved (k_inf = 0.4965). A forward and an adjoint problem are built on one cross-section object
    and solved; the object is then repopulated in place by LoadFromOpenSn with every cross section
    doubled, and both problems are solved again without being rebuilt. Their fluxes must be those
    of the new data, (2 (L - F))^-1 q and (2 (L - F)^T)^-1 q with L = diag(sigma_t) - S: (1, 0.6)
    and (0.6, 1) without fission. Also a one-group adjoint problem whose cross sections are
    repopulated by CreateSimpleOneGroup from sigma_t = 1 to sigma_t = 2, scattering ratio 0.5:
    phi = q, where the solve without scattering would give 0.5 q. Last, an adjoint problem on
    cross sections loaded by LoadFromOpenMC (HDPE.h5, 172 groups) and one on a three-group CEPXS
    file (LoadFromCEPXS) are solved, their cross sections are reloaded in place from the same file,
    and the solves are repeated without a rebuild: each flux must equal that of the same problem on
    a new object loaded from that file.
  tgdsa: the two-group medium of transport_tgdsa_two_group_infinite_medium.py, sigma_t = (1, 2),
    S = [[0.5, 0.09], [0.45, 1.9]], q = (1, 0), solved by classic Richardson with TGDSA alone (60
    iterations allowed) as there, in adjoint mode: an adjoint problem on its own cross sections, and
    an adjoint problem built first and then a forward problem on one shared object. The two-grid
    collapse and correction must use the transposed transfer matrix in the adjoint problems and
    the forward one in the forward problem; a forward collapse makes the adjoint iteration
    diverge. (WGDSA is left out: with WGDSA and TGDSA the adjoint iteration does not converge on
    this medium.)

Checks of the scenarios:
  SHARED_XS_KEIGEN_<FA|AF>_K_<FWD|ADJ>_REL_ERR: |k / k_inf - 1|.
  SHARED_XS_KEIGEN_<FA|AF>_RATIO_<FWD|ADJ>_REL_ERR: the relative error of P/R.
  SHARED_XS_CMFD_K_<FWD|ADJ>_REL_ERR, SHARED_XS_CMFD_RATIO_<FWD|ADJ>_REL_ERR: the same (FA order).
  SHARED_XS_CMFD_FWD_CONVERGED: 1 when the CMFD power iteration stopped before its iteration limit.
  SHARED_XS_CMFD_GROUPS_*: the SHARED_XS_CMFD_* checks of keigen_cmfd_groups.
  SHARED_XS_CMFD_ADJOINT_G<1|2>_<A|AF>_<K|RATIO>_REL_ERR and _CONVERGED: the same for the
    adjoint problem of keigen_cmfd_adjoint (G: group_aggregation_size).
  SHARED_XS_DELAYED_<FA|AF>_<K|RATIO>_FWD_REL_ERR: the k and P/R errors of keigen_delayed.
  SHARED_XS_FIRST_COLLISION_<FWD|ADJ>_REL_DIFF: max |phi - phi_unshared| / max |phi_unshared|.
  SHARED_XS_FIRST_COLLISION_FISSION_<FWD|ADJ>_REL_DIFF: the same for first_collision_fission.
  SHARED_XS_RELOAD_<FWD|ADJ>_ERR, SHARED_XS_RELOAD_FISSION_<FWD|ADJ>_ERR,
  SHARED_XS_RELOAD_ONE_GROUP_ADJ_ERR: max |phi / phi_exact - 1| after the reload.
  SHARED_XS_RELOAD_<OPENMC|CEPXS>_ADJ_REL_DIFF: max |phi - phi_new| / max |phi_new|, phi_new the
    flux on the new object.
  SHARED_XS_TGDSA_ADJ_ERR, SHARED_XS_TGDSA_AF_<ADJ|FWD>_ERR: max |phi / phi_exact - 1| of tgdsa.
"""

import os
import struct
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import PointSource, VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab, GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.solver import PowerIterationKEigenSolver, CMFDAcceleration
    from pyopensn.solver import UncollidedProblem, UncollidedSolver

SIGMA_T = np.array([1.0, 1.0])
SIGMA_A = np.array([0.2, 0.5])
# (from-group, to-group, value) of the isotropic transfer matrix
TRANSFERS = [(0, 0, 0.5), (0, 1, 0.3), (1, 1, 0.5)]

SIGMA_F = np.array([0.05, 0.25])
NU = np.array([2.5, 2.4])
CHI = np.array([0.9, 0.1])

# Delayed-neutron data of keigen_delayed (NU = NU_PROMPT + NU_DELAYED, CHI is the prompt spectrum)
NU_PROMPT = np.array([2.4, 2.3])
NU_DELAYED = np.array([0.1, 0.1])
CHI_DELAYED = [np.array([0.3, 0.7]), np.array([0.6, 0.4])]
DECAY_CONSTANTS = [0.08, 0.6]
FRACTIONAL_YIELDS = [0.4, 0.6]

# The medium of the tgdsa scenario
TGDSA_SIGMA_T = np.array([1.0, 2.0])
TGDSA_TRANSFER = np.array([[0.5, 0.09], [0.45, 1.9]])  # row: to-group

CMFD_MAX_ITERS = 2000

OPENMC_XS_FILE = "../../../../assets/xs/HDPE.h5"

BOUNDARIES = [
    {"name": "zmin", "type": "reflecting"},
    {"name": "zmax", "type": "reflecting"},
]


def scaled_sigma_f(scale, fission_scale):
    return SIGMA_F * (scale * fission_scale)


def section(name, values):
    rows = "".join(f"{g} {float(v)!r}\n" for g, v in enumerate(values))
    return f"{name}_BEGIN\n{rows}{name}_END\n"


def xs_text(fission=None, scale=1.0, fission_scale=1.0):
    """
    The medium in the OpenSn format, every cross section multiplied by scale and sigma_f also by
    fission_scale. fission is None, "prompt" or "delayed".
    """
    text = "NUM_GROUPS 2\nNUM_MOMENTS 1\n"
    if fission == "delayed":
        text += f"NUM_PRECURSORS {len(DECAY_CONSTANTS)}\n"
    text += section("SIGMA_T", scale * SIGMA_T) + section("SIGMA_A", scale * SIGMA_A)
    text += "TRANSFER_MOMENTS_BEGIN\n"
    for g_from, g_to, value in TRANSFERS:
        text += f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {scale * value!r}\n"
    text += "TRANSFER_MOMENTS_END\n"
    if fission is not None:
        text += section("SIGMA_F", scaled_sigma_f(scale, fission_scale))
    if fission == "prompt":
        text += section("NU", NU) + section("CHI", CHI)
    elif fission == "delayed":
        text += section("NU_PROMPT", NU_PROMPT) + section("NU_DELAYED", NU_DELAYED)
        text += section("CHI_PROMPT", CHI)
        text += "CHI_DELAYED_BEGIN\n"
        for j, spectrum in enumerate(CHI_DELAYED):
            for g, value in enumerate(spectrum):
                text += f"G_PRECURSOR_VAL {g} {j} {float(value)!r}\n"
        text += "CHI_DELAYED_END\n"
        text += section("PRECURSOR_DECAY_CONSTANTS", DECAY_CONSTANTS)
        text += section("PRECURSOR_FRACTIONAL_YIELDS", FRACTIONAL_YIELDS)
    return text


def tgdsa_xs_text():
    text = "NUM_GROUPS 2\nNUM_MOMENTS 1\n" + section("SIGMA_T", TGDSA_SIGMA_T)
    text += "TRANSFER_MOMENTS_BEGIN\n"
    for g_to in range(2):
        for g_from in range(2):
            text += f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {float(TGDSA_TRANSFER[g_to, g_from])!r}\n"
    return text + "TRANSFER_MOMENTS_END\n"


def load_xs(fission=None, scale=1.0, fission_scale=1.0, xs=None, text=None):
    """Load the medium (or text) into a new cross-section object, or into xs in place."""
    path = f"adjoint_mode_shared_xs_{scenario}_r{rank}.xs"
    with open(path, "w") as f:
        f.write(xs_text(fission, scale, fission_scale) if text is None else text)
    if xs is None:
        xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    os.remove(path)
    return xs


def source(q):
    return VolumetricSource(block_ids=[0], group_strength=list(q))


def groupset(quadrature, last_group=1):
    return {
        "groups_from_to": [0, last_group],
        "angular_quadrature": quadrature,
        "inner_linear_method": "petsc_gmres",
        "l_abs_tol": 1.0e-12,
        "l_max_its": 200,
    }


def options(adjoint):
    return {
        "adjoint": adjoint,
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": False,
    }


def slab_quadrature():
    return GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)


def build_problem(grid, xs, q, adjoint):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=len(q),
        groupsets=[groupset(slab_quadrature(), last_group=len(q) - 1)],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source(q)],
        boundary_conditions=BOUNDARIES,
        options=options(adjoint),
    )


def build_source_solver(problem):
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    return solver


def max_rel_err(problem, expected):
    phi = np.array(problem.GetPhiNewLocal()).reshape(-1, len(expected))
    local = float(np.max(np.abs(phi / np.array(expected) - 1.0))) if phi.size > 0 else 0.0
    return MPIAllReduce(local if np.isfinite(local) else 1.0e300, "max")


def fission_operators(fission="prompt", scale=1.0, fission_scale=1.0):
    """
    Return L = diag(sigma_t) - S and the production matrix F (row: to-group, delayed production
    included) of xs_text(fission, scale, fission_scale).
    """
    S = np.zeros((2, 2))
    for g_from, g_to, value in TRANSFERS:
        S[g_to, g_from] = scale * value
    L = np.diag(scale * SIGMA_T) - S
    sigma_f = scaled_sigma_f(scale, fission_scale)
    if fission is None:
        F = np.zeros((2, 2))
    elif fission == "prompt":
        F = np.outer(CHI, NU * sigma_f)
    else:
        chi_delayed = sum(y * chi for y, chi in zip(FRACTIONAL_YIELDS, CHI_DELAYED))
        F = np.outer(CHI, NU_PROMPT * sigma_f) + np.outer(chi_delayed, NU_DELAYED * sigma_f)
    return L, F


def dominant(L, F):
    """Dominant eigenvalue and eigenvector of L^-1 F."""
    values, vectors = np.linalg.eig(np.linalg.solve(L, F))
    i = np.argmax(values.real)
    return float(values[i].real), vectors[:, i].real


def production_ratio_reference(adjoint, fission="prompt"):
    L, F = fission_operators(fission)
    if adjoint:
        L, F = L.T, F.T
    _, phi = dominant(L, F)
    return float(np.sum(F @ phi) / (SIGMA_F @ phi))


def build_keigen(grid, xs, adjoint, cmfd=False, group_aggregation_size=1):
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[groupset(slab_quadrature())],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=BOUNDARIES,
        options=options(adjoint),
    )
    if cmfd:
        # The settings of keigenvalue_transport_1d_1g_cmfd.py.
        acceleration = CMFDAcceleration(
            problem=problem,
            current_closure="auto",
            coarse_mesh="local_aggregation",
            aggregation_size=2,
            group_aggregation_size=group_aggregation_size,
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
            problem=problem, acceleration=acceleration, max_iters=CMFD_MAX_ITERS, k_tol=1.0e-10
        )
    else:
        solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-12)
    solver.Initialize()
    return problem, solver


def keigen_errors(problem, solver, adjoint, fission="prompt"):
    solver.Execute()
    k_inf = dominant(*fission_operators(fission))[0]
    ratio = problem.ComputeFissionProduction("new") / problem.ComputeFissionRate("new")
    return (
        abs(solver.GetEigenvalue() / k_inf - 1.0),
        abs(ratio / production_ratio_reference(adjoint, fission) - 1.0),
    )


def require_procs(num_procs):
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")


def slab_grid():
    nodes = [0.1 * i for i in range(11)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    return grid


def run_fixed_source(results):
    require_procs(2)
    grid = slab_grid()

    # A. One cross-section object shared by a forward and an adjoint problem.
    xs_shared = load_xs()
    forward = build_problem(grid, xs_shared, (1.0, 0.0), adjoint=False)
    forward_solver = build_source_solver(forward)
    forward_solver.Execute()
    results["SHARED_XS_FWD_FIRST_ERR"] = max_rel_err(forward, (2.0, 1.2))

    adjoint = build_problem(grid, xs_shared, (0.0, 1.0), adjoint=True)
    adjoint_solver = build_source_solver(adjoint)
    adjoint_solver.Execute()
    results["SHARED_XS_ADJ_ERR"] = max_rel_err(adjoint, (1.2, 2.0))

    forward_solver.Execute()
    results["SHARED_XS_FWD_AGAIN_ERR"] = max_rel_err(forward, (2.0, 1.2))
    adjoint_solver.Execute()
    results["SHARED_XS_ADJ_AGAIN_ERR"] = max_rel_err(adjoint, (1.2, 2.0))

    # B. Cross sections scaled between two switches into adjoint mode.
    xs_scaled = load_xs()
    phys = build_problem(grid, xs_scaled, (1.0, 0.0), adjoint=False)
    solver = build_source_solver(phys)
    phys.SetAdjoint(True)
    phys.SetAdjoint(False)
    xs_scaled.Scale(2.0)
    phys.SetAdjoint(True)
    phys.SetBoundaryOptions(boundary_conditions=BOUNDARIES)
    phys.SetVolumetricSources(volumetric_sources=[source((0.0, 1.0))])
    solver.Execute()
    results["SCALED_XS_ADJ_ERR"] = max_rel_err(phys, (0.6, 1.0))

    # C. Cross sections scaled under a live adjoint problem after a forward problem is built.
    xs_live = load_xs()
    live_adjoint = build_problem(grid, xs_live, (0.0, 1.0), adjoint=True)
    live_adjoint_solver = build_source_solver(live_adjoint)
    live_forward = build_problem(grid, xs_live, (1.0, 0.0), adjoint=False)
    live_forward_solver = build_source_solver(live_forward)
    xs_live.Scale(2.0)
    live_adjoint_solver.Execute()
    results["SHARED_XS_LIVE_ADJ_SCALED_ERR"] = max_rel_err(live_adjoint, (0.6, 1.0))
    live_forward_solver.Execute()
    results["SHARED_XS_LIVE_FWD_SCALED_ERR"] = max_rel_err(live_forward, (1.0, 0.6))


def write_cepxs_file():
    """Write a three-group CEPXS file (the fixture of csda_balance.py without the threshold)."""
    path = f"adjoint_mode_shared_xs_{scenario}_r{rank}.bxslib"
    with open(path, "wb") as f:

        def record(payload):
            marker = struct.pack("<I", len(payload))
            f.write(marker + payload + marker)

        scattering = [[0.20, 0.0, 0.0], [0.10, 0.15, 0.0], [0.03, 0.12, 0.10]]
        table = [0.0] * 33
        for g in range(3):
            table[11 * g + 1] = [0.05, 0.03, 0.01][g]
            table[11 * g + 2] = [0.5, 0.35, 0.2][g]
            table[11 * g + 4] = [0.4, 0.3, 0.2][g]
            table[11 * g + 7] = [1.0, 1.2, 1.5][g]
            for gp in range(g + 1):
                table[11 * g + 8 + g - gp] = scattering[g][gp]
        record(b"ADJOINT MODE SHARED XS")
        record(struct.pack("<8i", 3, 1, 11, 8, 9, 1, 0, 1))
        record(struct.pack("<4d", 3.0, 2.0, 1.0, 0.0))
        record(struct.pack("<33d", *table))
    return path


def reload_rel_diff(grid, load):
    """
    Solve an adjoint problem on load(xs), reload its cross sections in place with load(xs), solve
    again without a rebuild, and compare with the same problem on load(MultiGroupXS()).
    """

    def solve(xs):
        q = [1.0] * xs.num_groups
        problem = build_problem(grid, xs, q, adjoint=True)
        solver = build_source_solver(problem)
        solver.Execute()
        return problem, solver

    xs = MultiGroupXS()
    load(xs)
    problem, solver = solve(xs)
    load(xs)
    solver.Execute()
    phi = np.array(problem.GetPhiNewLocal())
    new_problem, _ = solve(load(MultiGroupXS()))
    phi_new = np.array(new_problem.GetPhiNewLocal())
    local_diff = float(np.max(np.abs(phi - phi_new))) if phi.size > 0 else 0.0
    local_max = float(np.max(np.abs(phi_new))) if phi.size > 0 else 0.0
    diff = MPIAllReduce(local_diff if np.isfinite(local_diff) else 1.0e300, "max")
    return diff / MPIAllReduce(local_max, "max")


def run_reload(results):
    require_procs(2)
    grid = slab_grid()
    q_forward, q_adjoint = np.array([1.0, 0.0]), np.array([0.0, 1.0])

    for fission, prefix in ((None, "SHARED_XS_RELOAD"), ("prompt", "SHARED_XS_RELOAD_FISSION")):
        xs = load_xs(fission, fission_scale=0.5)
        forward = build_problem(grid, xs, q_forward, adjoint=False)
        forward_solver = build_source_solver(forward)
        forward_solver.Execute()
        adjoint = build_problem(grid, xs, q_adjoint, adjoint=True)
        adjoint_solver = build_source_solver(adjoint)
        adjoint_solver.Execute()

        load_xs(fission, scale=2.0, fission_scale=0.5, xs=xs)
        L, F = fission_operators(fission, scale=2.0, fission_scale=0.5)
        adjoint_solver.Execute()
        results[f"{prefix}_ADJ_ERR"] = max_rel_err(adjoint, np.linalg.solve((L - F).T, q_adjoint))
        forward_solver.Execute()
        results[f"{prefix}_FWD_ERR"] = max_rel_err(forward, np.linalg.solve(L - F, q_forward))

    xs_one_group = MultiGroupXS()
    xs_one_group.CreateSimpleOneGroup(1.0, 0.5)
    one_group = build_problem(grid, xs_one_group, (1.0,), adjoint=True)
    one_group_solver = build_source_solver(one_group)
    one_group_solver.Execute()
    xs_one_group.CreateSimpleOneGroup(2.0, 0.5)
    one_group_solver.Execute()
    results["SHARED_XS_RELOAD_ONE_GROUP_ADJ_ERR"] = max_rel_err(one_group, (1.0,))

    def load_openmc(xs):
        xs.LoadFromOpenMC(OPENMC_XS_FILE, "set1", 294.0)
        return xs

    cepxs_path = write_cepxs_file()

    def load_cepxs(xs):
        xs.LoadFromCEPXS(cepxs_path, material_id=0, csda_format=False)
        return xs

    results["SHARED_XS_RELOAD_OPENMC_ADJ_REL_DIFF"] = reload_rel_diff(grid, load_openmc)
    results["SHARED_XS_RELOAD_CEPXS_ADJ_REL_DIFF"] = reload_rel_diff(grid, load_cepxs)
    os.remove(cepxs_path)


def run_tgdsa(results):
    require_procs(2)
    nodes = [1.0 * i for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    q = np.array([1.0, 0.0])
    L = np.diag(TGDSA_SIGMA_T) - TGDSA_TRANSFER

    def build(xs, adjoint):
        problem = DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=2,
            groupsets=[
                {
                    "groups_from_to": [0, 1],
                    "angular_quadrature": GLProductQuadrature1DSlab(
                        n_polar=16, scattering_order=0
                    ),
                    "inner_linear_method": "classic_richardson",
                    "l_abs_tol": 1.0e-10,
                    "l_max_its": 60,
                    "apply_wgdsa": False,
                    "apply_tgdsa": True,
                    "tgdsa_l_abs_tol": 1.0e-12,
                }
            ],
            xs_map=[{"block_ids": [0], "xs": xs}],
            volumetric_sources=[source(q)],
            boundary_conditions=BOUNDARIES,
            options=options(adjoint),
        )
        return problem, build_source_solver(problem)

    adjoint, adjoint_solver = build(load_xs(text=tgdsa_xs_text()), adjoint=True)
    adjoint_solver.Execute()
    results["SHARED_XS_TGDSA_ADJ_ERR"] = max_rel_err(adjoint, np.linalg.solve(L.T, q))

    xs = load_xs(text=tgdsa_xs_text())
    shared_adjoint, shared_adjoint_solver = build(xs, adjoint=True)
    shared_forward, shared_forward_solver = build(xs, adjoint=False)
    shared_forward_solver.Execute()
    results["SHARED_XS_TGDSA_AF_FWD_ERR"] = max_rel_err(shared_forward, np.linalg.solve(L, q))
    shared_adjoint_solver.Execute()
    results["SHARED_XS_TGDSA_AF_ADJ_ERR"] = max_rel_err(shared_adjoint, np.linalg.solve(L.T, q))


def run_keigen(results):
    require_procs(2)
    grid = slab_grid()

    xs_fa = load_xs("prompt")
    forward = build_keigen(grid, xs_fa, adjoint=False)
    adjoint = build_keigen(grid, xs_fa, adjoint=True)
    k_err, ratio_err = keigen_errors(*forward, adjoint=False)
    results["SHARED_XS_KEIGEN_FA_K_FWD_REL_ERR"] = k_err
    results["SHARED_XS_KEIGEN_FA_RATIO_FWD_REL_ERR"] = ratio_err
    k_err, ratio_err = keigen_errors(*adjoint, adjoint=True)
    results["SHARED_XS_KEIGEN_FA_K_ADJ_REL_ERR"] = k_err
    results["SHARED_XS_KEIGEN_FA_RATIO_ADJ_REL_ERR"] = ratio_err

    xs_af = load_xs("prompt")
    adjoint = build_keigen(grid, xs_af, adjoint=True)
    forward = build_keigen(grid, xs_af, adjoint=False)
    k_err, ratio_err = keigen_errors(*adjoint, adjoint=True)
    results["SHARED_XS_KEIGEN_AF_K_ADJ_REL_ERR"] = k_err
    results["SHARED_XS_KEIGEN_AF_RATIO_ADJ_REL_ERR"] = ratio_err
    k_err, ratio_err = keigen_errors(*forward, adjoint=False)
    results["SHARED_XS_KEIGEN_AF_K_FWD_REL_ERR"] = k_err
    results["SHARED_XS_KEIGEN_AF_RATIO_FWD_REL_ERR"] = ratio_err


def run_keigen_cmfd(results, group_aggregation_size=1, prefix="SHARED_XS_CMFD"):
    require_procs(2)
    grid = slab_grid()

    xs = load_xs("prompt")
    forward = build_keigen(
        grid, xs, adjoint=False, cmfd=True, group_aggregation_size=group_aggregation_size
    )
    adjoint = build_keigen(grid, xs, adjoint=True)
    k_err, ratio_err = keigen_errors(*forward, adjoint=False)
    results[f"{prefix}_K_FWD_REL_ERR"] = k_err
    results[f"{prefix}_RATIO_FWD_REL_ERR"] = ratio_err
    iterations = forward[1].GetNumPowerIterations()
    results[f"{prefix}_FWD_POWER_ITERATIONS"] = iterations
    results[f"{prefix}_FWD_CONVERGED"] = 1.0 if iterations < CMFD_MAX_ITERS else 0.0
    k_err, ratio_err = keigen_errors(*adjoint, adjoint=True)
    results[f"{prefix}_K_ADJ_REL_ERR"] = k_err
    results[f"{prefix}_RATIO_ADJ_REL_ERR"] = ratio_err


def run_keigen_cmfd_adjoint(results):
    require_procs(2)
    grid = slab_grid()

    for group_aggregation_size in (1, 2):
        for order in ("A", "AF"):
            xs = load_xs("prompt")
            problems = [
                build_keigen(
                    grid, xs, adjoint=True, cmfd=True, group_aggregation_size=group_aggregation_size
                )
            ]
            if order == "AF":
                problems.append(build_keigen(grid, xs, adjoint=False))
            prefix = f"SHARED_XS_CMFD_ADJOINT_G{group_aggregation_size}_{order}"
            k_err, ratio_err = keigen_errors(*problems[0], adjoint=True)
            results[f"{prefix}_K_REL_ERR"] = k_err
            results[f"{prefix}_RATIO_REL_ERR"] = ratio_err
            iterations = problems[0][1].GetNumPowerIterations()
            results[f"{prefix}_POWER_ITERATIONS"] = iterations
            results[f"{prefix}_CONVERGED"] = 1.0 if iterations < CMFD_MAX_ITERS else 0.0


def run_keigen_delayed(results):
    require_procs(2)
    grid = slab_grid()

    for order, modes in (("FA", (False, True)), ("AF", (True, False))):
        xs = load_xs("delayed")
        problems = {adjoint: build_keigen(grid, xs, adjoint=adjoint) for adjoint in modes}
        k_err, ratio_err = keigen_errors(*problems[False], adjoint=False, fission="delayed")
        results[f"SHARED_XS_DELAYED_{order}_K_FWD_REL_ERR"] = k_err
        results[f"SHARED_XS_DELAYED_{order}_RATIO_FWD_REL_ERR"] = ratio_err


def run_first_collision(results, fission=None, prefix="SHARED_XS_FIRST_COLLISION"):
    require_procs(1)
    nodes = [-1.0, -0.5, 0.0, 0.5, 1.0]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    quadrature = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=8, scattering_order=0)

    xs = load_xs(fission)
    file_name = f"adjoint_mode_shared_xs_{scenario}.h5"
    if os.path.exists(file_name):
        os.remove(file_name)
    uncollided = UncollidedProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[{"groups_from_to": [0, 1]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=[0.14, -0.19, 0.0], strength=[1.0, 0.5])],
        near_source=[RPPLogicalVolume(xmin=-1.01, xmax=1.01, ymin=-1.01, ymax=1.01, infz=True)],
        scattering_order=0,
    )
    uncollided_solver = UncollidedSolver(problem=uncollided, file_name=file_name)
    uncollided_solver.Initialize()
    uncollided_solver.Execute()

    def build(cross_sections, adjoint):
        if adjoint:
            settings = {"volumetric_sources": [source((0.0, 1.0))]}
        else:
            settings = {"uncollided_flux": file_name}
        problem = DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=2,
            groupsets=[groupset(quadrature)],
            xs_map=[{"block_ids": [0], "xs": cross_sections}],
            options=options(adjoint),
            **settings,
        )
        return problem, build_source_solver(problem)

    def solve(problem, solver):
        solver.Execute()
        return np.array(problem.GetPhiNewLocal())

    def rel_diff(phi, phi_reference):
        return float(np.max(np.abs(phi - phi_reference)) / np.max(np.abs(phi_reference)))

    adjoint = build(xs, adjoint=True)
    forward = build(xs, adjoint=False)
    phi_forward = solve(*forward)
    phi_adjoint = solve(*adjoint)
    phi_forward_unshared = solve(*build(load_xs(fission), adjoint=False))
    phi_adjoint_unshared = solve(*build(load_xs(fission), adjoint=True))
    os.remove(file_name)
    results[f"{prefix}_FWD_REL_DIFF"] = rel_diff(phi_forward, phi_forward_unshared)
    results[f"{prefix}_ADJ_REL_DIFF"] = rel_diff(phi_adjoint, phi_adjoint_unshared)


if __name__ == "__main__":
    if "scenario" not in globals():
        scenario = "fixed_source"

    scenarios = {
        "fixed_source": run_fixed_source,
        "reload": run_reload,
        "tgdsa": run_tgdsa,
        "keigen": run_keigen,
        "keigen_cmfd": run_keigen_cmfd,
        "keigen_cmfd_groups": lambda r: run_keigen_cmfd(r, 2, "SHARED_XS_CMFD_GROUPS"),
        "keigen_cmfd_adjoint": run_keigen_cmfd_adjoint,
        "keigen_delayed": run_keigen_delayed,
        "first_collision": run_first_collision,
        "first_collision_fission": lambda r: run_first_collision(
            r, "prompt", "SHARED_XS_FIRST_COLLISION_FISSION"
        ),
    }
    results = {}
    scenarios[scenario](results)

    if rank == 0:
        for key, value in results.items():
            print(f"{key}={value:.6e}")
