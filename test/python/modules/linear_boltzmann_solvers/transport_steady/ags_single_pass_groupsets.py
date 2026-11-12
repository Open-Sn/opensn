#!/usr/bin/env python3
"""
Across-groupset (AGS) iteration with and without coupling to earlier groupsets.

A 1D slab with the 168-group library and a source in the first group is solved with all groups
in one groupset, the reference, and split into two groupsets in two ways. The library scatters up
from groups 70 and above into groups 69 and above.

- Split A puts all of the upscattering inside the second groupset (0-62, 63-167), so no groupset
  receives a source from one solved after it, and one pass over the groupsets solves the
  problem: with max_ags_iterations = 1 the solve converges.
- Split B splits the upscattering range (0-99, 100-167), so the second groupset feeds the first,
  and one pass is not enough: with max_ags_iterations = 1 the solve does not converge. With the
  default limit and a tight AGS tolerance it converges to the reference.

Both splits are compared with the reference group by group, with inner tolerances of 1e-12. The
difference is largest, about 2e-9, in the thermal groups with 1e-3 of the peak flux, where the
slow convergence of the upscattering amplifies the inner error; it falls in proportion to the
inner tolerance. Split B is also solved with the source scaled by 1e-12, with an inner tolerance
of 1e-6 and the relative (pointwise) AGS check: the problem is linear, so its flux must be the
scaled flux of the unscaled solve, iteration for iteration, up to round-off that the Krylov
iterations amplify to about 1e-8. A groupset solve that discarded a small initial guess would
restart and change the flux at the level of the inner tolerance, about 1e-5.
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
    from pyopensn.post import VolumePostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

NUM_GROUPS = 168


def solve(splits, scale=1.0, l_abs_tol=1.0e-12, options=None):
    """Return the group fluxes integrated over the slab and whether the solve converged."""
    nodes = [30.0 * i / 60 for i in range(61)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/xs_168g.xs")
    quadrature = GLProductQuadrature1DSlab(n_polar=16, scattering_order=1)
    groupsets = [{"groups_from_to": split, "angular_quadrature": quadrature,
                  "inner_linear_method": "petsc_gmres", "l_abs_tol": l_abs_tol,
                  "l_max_its": 500, "gmres_restart_interval": 100} for split in splits]
    source = VolumetricSource(block_ids=[0],
                              group_strength=[scale] + [0.0] * (NUM_GROUPS - 1))
    problem = DiscreteOrdinatesProblem(mesh=grid, num_groups=NUM_GROUPS, groupsets=groupsets,
                                       xs_map=[{"block_ids": [0], "xs": xs}],
                                       volumetric_sources=[source], options=options or {})
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    post = VolumePostprocessor(problem=problem, value_type="integral")
    post.Execute()
    flux = np.asarray(post.GetValue()).reshape(-1, NUM_GROUPS).sum(axis=0)
    return flux, solver.IsConverged()


def max_relative_difference(a, b):
    """Largest difference relative to b, over the groups; groups where b is zero compare the
    absolute difference, which must be zero too."""
    return float(np.max(np.abs(a - b) / np.where(b != 0.0, np.abs(b), 1.0)))


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    split_a = [(0, 62), (63, NUM_GROUPS - 1)]
    split_b = [(0, 99), (100, NUM_GROUPS - 1)]
    reference, _ = solve([(0, NUM_GROUPS - 1)])
    a_one_pass, a_converged = solve(split_a, options={"max_ags_iterations": 1})
    _, b_one_pass_converged = solve(split_b, options={"max_ags_iterations": 1})
    pointwise = {"ags_convergence_check": "pointwise", "ags_tolerance": 1.0e-12}
    b_converged_flux, b_converged = solve(split_b, options=pointwise)
    loose = {"ags_convergence_check": "pointwise", "ags_tolerance": 1.0e-6}
    b_unit, _ = solve(split_b, l_abs_tol=1.0e-6, options=loose)
    b_scaled, _ = solve(split_b, scale=1.0e-12, l_abs_tol=1.0e-6, options=loose)

    if rank == 0:
        print(f"AGS_SPLIT_A_ONE_PASS_CONVERGED={int(a_converged)}")
        print(f"AGS_SPLIT_A_MAX_REL_DIFF={max_relative_difference(a_one_pass, reference):.6e}")
        print(f"AGS_SPLIT_B_ONE_PASS_CONVERGED={int(b_one_pass_converged)}")
        print(f"AGS_SPLIT_B_CONVERGED={int(b_converged)}")
        print(f"AGS_SPLIT_B_MAX_REL_DIFF="
              f"{max_relative_difference(b_converged_flux, reference):.6e}")
        print(f"AGS_SOURCE_SCALE_MAX_REL_DIFF="
              f"{max_relative_difference(b_scaled / 1.0e-12, b_unit):.6e}")
