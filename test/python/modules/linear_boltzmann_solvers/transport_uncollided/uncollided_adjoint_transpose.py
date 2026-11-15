#!/usr/bin/env python3
"""Verify that an adjoint first-collision solve is the forward solve with transposed cross sections.

The adjoint transport operator with transfer matrices S and production matrix F is the forward
operator with S^T and F^T, with the directions reversed. OpenSn solves the adjoint along the
quadrature directions with the transposed matrices and then reorients the solution, which flips
the sign of the flux moments of odd degree. The uncollided flux of a point source depends only on
sigma_t, which is the same in both modes. So the adjoint flux moments of a problem with S, F, and
an adjoint uncollided file from a point detector must equal (-1)^ell times the flux moments of a
forward problem with S^T, F^T, and a forward uncollided file from a point source with the same
location and strengths.

The cross sections have two groups, P1 scattering with up- and downscatter, and fission, with
nonsymmetric S and F. The case checks that the adjoint first-collision source uses S^T and F^T,
and that the uncollided moments are added before the adjoint solution is reoriented. Adding them
after the reorientation gives P1 moments of the wrong sign near the detector. The independent
reference for the adjoint uncollided path is the single-scatter integral of
uncollided_adjoint_reciprocity.py.

Checks:
  UncollidedAdjointTransposeMaxRelDiff: max |(-1)^ell phi_adj - phi_fwd^T| / max |phi_fwd^T| (0 to
    round-off).
  UncollidedAdjointTransposeSensitive: 1 if the forward solution with S and F differs from the one
    with S^T and F^T by more than 1e-2 relative, and the P1 moments are larger than 1e-2 of the
    scalar flux, so that the comparison can detect the untransposed matrices and the wrong sign.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier

    def global_max(value):
        return MPI.COMM_WORLD.allreduce(value, op=MPI.MAX)

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
        UncollidedProblem,
        UncollidedSolver,
    )
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier

    def global_max(value):
        return MPIAllReduce(value, "max")


SIGMA_T = (1.0, 0.8)
# S[ell][to][from], nonsymmetric, with up- and downscatter.
S = (((0.30, 0.05), (0.25, 0.40)),
     ((0.10, 0.00), (0.08, 0.12)))
# F[to][from], nonsymmetric. The reader derives nu from the column sums of the production matrix and
# sigma_f, and both F and F^T give nu > 1 in every group. The infinite-medium k is 0.38.
F = ((0.06, 0.20), (0.05, 0.0))
SIGMA_F = (0.02, 0.04)
POINT = (1.27, 1.61, 1.43)
STRENGTH = [0.3, 1.0]
NUM_GROUPS = 2
NUM_MOMENTS = 4


def remove(file_name):
    if rank == 0 and os.path.exists(file_name):
        os.remove(file_name)
    barrier()


def write_xs(file_name, transpose):
    if rank == 0:
        lines = ["NUM_GROUPS 2", "NUM_MOMENTS 2", "SIGMA_T_BEGIN"]
        lines += [f"{g} {SIGMA_T[g]}" for g in range(2)]
        lines += ["SIGMA_T_END", "SIGMA_F_BEGIN"]
        lines += [f"{g} {SIGMA_F[g]}" for g in range(2)]
        lines += ["SIGMA_F_END", "CHI_BEGIN", "0 1.0", "1 0.0", "CHI_END"]
        lines += ["PRODUCTION_MATRIX_BEGIN"]
        for g_to in range(2):
            for g_from in range(2):
                value = F[g_from][g_to] if transpose else F[g_to][g_from]
                lines.append(f"GPRIME_G_VAL {g_from} {g_to} {value}")
        lines += ["PRODUCTION_MATRIX_END", "TRANSFER_MOMENTS_BEGIN"]
        for ell in range(2):
            for g_to in range(2):
                for g_from in range(2):
                    value = S[ell][g_from][g_to] if transpose else S[ell][g_to][g_from]
                    lines.append(f"M_GFROM_GTO_VAL {ell} {g_from} {g_to} {value}")
        lines += ["TRANSFER_MOMENTS_END"]
        with open(file_name, "w") as f:
            f.write("\n".join(lines) + "\n")
    barrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(file_name)
    remove(file_name)
    return xs


def solve(grid, xs, adjoint, tag):
    file_name = f"uncollided_adjoint_transpose_{tag}.h5"
    remove(file_name)
    generator = UncollidedProblem(
        mesh=grid,
        num_groups=NUM_GROUPS,
        groupsets=[{"groups_from_to": [0, 1]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=list(POINT), strength=STRENGTH)],
        near_source=[
            RPPLogicalVolume(
                xmin=POINT[0] - 0.6,
                xmax=POINT[0] + 0.6,
                ymin=POINT[1] - 0.6,
                ymax=POINT[1] + 0.6,
                zmin=POINT[2] - 0.6,
                zmax=POINT[2] + 0.6,
            )
        ],
        scattering_order=1,
        options={"adjoint": adjoint},
    )
    generation = UncollidedSolver(problem=generator, file_name=file_name)
    generation.Initialize()
    generation.Execute()

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=NUM_GROUPS,
        groupsets=[
            {
                "groups_from_to": [0, 1],
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=4, n_azimuthal=8, scattering_order=1
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        uncollided_flux=file_name,
        options={"adjoint": adjoint, "verbose_inner_iterations": False},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    remove(file_name)
    return list(problem.GetPhiNewLocal())


def moment_sign(k):
    # Local flux moments are stored by node, then moment, then group. The 3D P1 moments are (0, 0),
    # (1, -1), (1, 0), and (1, 1).
    return 1.0 if (k // NUM_GROUPS) % NUM_MOMENTS == 0 else -1.0


def is_p1(k):
    return (k // NUM_GROUPS) % NUM_MOMENTS != 0


if __name__ == "__main__":
    nodes = [0.5 * i for i in range(7)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs_s = write_xs("uncollided_adjoint_transpose_S.xs", transpose=False)
    xs_st = write_xs("uncollided_adjoint_transpose_ST.xs", transpose=True)

    phi_adj = solve(grid, xs_s, adjoint=True, tag="adjoint")
    phi_fwd_t = solve(grid, xs_st, adjoint=False, tag="forward_transposed")
    phi_fwd = solve(grid, xs_s, adjoint=False, tag="forward")

    scale = global_max(max(abs(v) for v in phi_fwd_t))
    diff = global_max(
        max(abs(moment_sign(k) * a - b) for k, (a, b) in enumerate(zip(phi_adj, phi_fwd_t)))
    )
    transpose_effect = global_max(max(abs(a - b) for a, b in zip(phi_fwd, phi_fwd_t)))
    scalar_scale = global_max(max(abs(v) for k, v in enumerate(phi_fwd_t) if not is_p1(k)))
    p1_scale = global_max(max(abs(v) for k, v in enumerate(phi_fwd_t) if is_p1(k)))
    sensitive = transpose_effect > 1.0e-2 * scale and p1_scale > 1.0e-2 * scalar_scale

    if rank == 0:
        print(f"UncollidedAdjointTransposeTransposeEffect={transpose_effect / scale:.8e}")
        print(f"UncollidedAdjointTransposeP1Ratio={p1_scale / scalar_scale:.8e}")
        print(f"UncollidedAdjointTransposeMaxRelDiff={diff / scale:.8e}")
        print(f"UncollidedAdjointTransposeSensitive={int(sensitive)}")
        sys.stdout.flush()
