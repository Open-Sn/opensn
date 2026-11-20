#!/usr/bin/env python3
"""
GetPhiNewLocal and GetPhiOldLocal return flux moments in node, moment, group order for any
groupset split.

A 1D slab with reflecting boundaries on both ends and a uniform source is an infinite medium: the
flux is flat, the first angular moment is zero, and the scalar flux solves the group balance

    sigma_t phi_g = sum_g' sigma_s(g' -> g) phi_g' + Q_g.

The library simple_upscatter.xs has sigma_t = 1 in three groups, a within-group transfer of 0.05,
and a transfer of 0.45 to each other group. The balance matrix is then 1.4 I - 0.45 J, with J the
matrix of ones, whose inverse is (I + 9 J) / 1.4, so

    phi_g = (Q_g + 9 sum(Q)) / 1.4.

With Q = (1, 2, 0.5) this is (23.2142857, 23.9285714, 22.8571429). The test evaluates it with
numpy from the transfer matrix above, after checking that matrix against the library's sigma_t and
sigma_a. The group fluxes differ by about 3%, so a view in any other order fails the comparison.

The problem is solved with all groups in one groupset, where the solver stores flux moments in
node, moment, group order, and with two groupsets of unequal sizes, (0) and (1, 2), where it stores
each groupset separately and the accessors return a reordered copy. The upscatter couples the two
groupsets, so the second case also runs across-groupset iterations. The inner and across-groupset
tolerances are 1e-12; with a scattering ratio of 0.95 the iteration error is below about 1e-10,
so the accepted relative error is 1e-8. The first moment is checked against 1e-8 in absolute value.

Both views must be read-only. After the two-groupset problem is solved again with the source
doubled, a view taken before must keep the old values (it is a copy) and a new view must hold twice
the old ones (the problem is linear).

Each rank checks its own local vectors and prints its own results.
"""

import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

NUM_GROUPS = 3
NUM_MOMENTS = 2
SOURCE = [1.0, 2.0, 0.5]


# Zeroth transfer moment of simple_upscatter.xs, TRANSFER[g_from][g_to].
TRANSFER = np.array([[0.05, 0.45, 0.45], [0.45, 0.05, 0.45], [0.45, 0.45, 0.05]])


def analytic_flux(xs):
    """Infinite-medium scalar flux: solves (diag(sigma_t) - TRANSFER^T) phi = Q."""
    sigma_t = np.asarray(xs.sigma_t)
    # The library's absorption is its total minus its outscatter, which checks TRANSFER.
    if not np.allclose(np.asarray(xs.sigma_a), sigma_t - TRANSFER.sum(axis=1), rtol=0, atol=1e-14):
        raise RuntimeError("TRANSFER does not match simple_upscatter.xs.")
    return np.linalg.solve(np.diag(sigma_t) - TRANSFER.T, np.asarray(SOURCE))


def make_problem(xs, splits):
    nodes = [4.0 * i / 8 for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    quadrature = GLProductQuadrature1DSlab(n_polar=4, scattering_order=1)
    groupsets = [{"groups_from_to": split, "angular_quadrature": quadrature,
                  "inner_linear_method": "petsc_gmres", "l_abs_tol": 1.0e-12, "l_max_its": 300}
                 for split in splits]
    return DiscreteOrdinatesProblem(
        mesh=grid, num_groups=NUM_GROUPS, groupsets=groupsets,
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=SOURCE)],
        boundary_conditions=[{"name": "zmin", "type": "reflecting"},
                             {"name": "zmax", "type": "reflecting"}],
        options={"max_ags_iterations": 500, "ags_tolerance": 1.0e-12,
                 "ags_convergence_check": "pointwise",
                 "verbose_inner_iterations": False, "verbose_outer_iterations": False})


def solve(problem):
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()


def moments(view):
    """Reshape a flux view to [node, moment, group], the documented order."""
    return np.asarray(view).reshape(-1, NUM_MOMENTS, NUM_GROUPS)


def order_errors(view, phi_exact):
    """Largest relative error of the scalar flux and largest first moment."""
    values = moments(view)
    phi0_err = float(np.max(np.abs(values[:, 0, :] / phi_exact - 1.0))) if values.size else 0.0
    phi1_max = float(np.max(np.abs(values[:, 1, :]))) if values.size else 0.0
    return phi0_err, phi1_max


def is_read_only(view):
    try:
        view[0] = 0.0
    except TypeError:
        return True
    return False


if __name__ == "__main__":
    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/simple_upscatter.xs")
    phi_exact = analytic_flux(xs)

    phi0_err = 0.0
    phi1_max = 0.0
    read_only = True

    single = make_problem(xs, [(0, 2)])
    solve(single)
    for view in (single.GetPhiNewLocal(), single.GetPhiOldLocal()):
        e0, e1 = order_errors(view, phi_exact)
        phi0_err, phi1_max = max(phi0_err, e0), max(phi1_max, e1)
        read_only = read_only and (len(view) == 0 or is_read_only(view))

    split = make_problem(xs, [(0, 0), (1, 2)])
    solve(split)
    before = split.GetPhiNewLocal()
    for view in (before, split.GetPhiOldLocal()):
        e0, e1 = order_errors(view, phi_exact)
        phi0_err, phi1_max = max(phi0_err, e0), max(phi1_max, e1)
        read_only = read_only and (len(view) == 0 or is_read_only(view))
    before_values = np.array(before, copy=True)

    split.SetVolumetricSources(
        volumetric_sources=[VolumetricSource(block_ids=[0],
                                             group_strength=[2.0 * q for q in SOURCE])],
        clear_volumetric_sources=True)
    solve(split)
    after = split.GetPhiNewLocal()
    e0, e1 = order_errors(after, 2.0 * phi_exact)
    phi0_err, phi1_max = max(phi0_err, e0), max(phi1_max, e1)
    snapshot = bool(np.array_equal(np.asarray(before), before_values))

    print(f"PHI_LOCAL_ORDER_RANK{rank}_PHI0_MAX_REL_ERR={phi0_err:.6e}")
    print(f"PHI_LOCAL_ORDER_RANK{rank}_PHI1_MAX_ABS={phi1_max:.6e}")
    print(f"PHI_LOCAL_ORDER_RANK{rank}_READ_ONLY={int(read_only)}")
    print(f"PHI_LOCAL_ORDER_RANK{rank}_SPLIT_VIEW_IS_SNAPSHOT={int(snapshot)}")
