#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
WGDSA and TGDSA iterates of an adjoint problem equal those of the forward problem with transposed
scattering (2 ranks).

An adjoint problem with transfer matrices S sweeps along the quadrature directions with S^T, so
every iterate must equal the iterate of a forward problem given S^T directly: the same sweeps, the
same WGDSA correction with the residual of the lagged angular fluxes (an opposing reflecting
boundary, and cycles broken within and between partitions), and the same TGDSA spectrum. The
diffusion coefficients enter every correction. They use the transport correction of the P1
scattering, which is the in-scatter sum_g' S1[g][g'] for the adjoint and the out-scatter for the
forward problem with S^T, the same values. The scattering is P1 with nonsymmetric up- and
downscatter, so the forward coefficients of S differ.

Each case runs exactly 5 inner iterations, so the comparison is between unconverged accelerated
iterates. The adjoint flux moments are in the physical orientation, in which the P1 moments change
sign, so the adjoint moments of degree ell are multiplied by (-1)^ell before the comparison.
  inf_xy:        10x10 square, all faces reflecting (opposing pairs lag), uniform source.
  cut_row_mixed: 10x9 square, LinearGraphPartitioner (on 2 ranks the cut has a step, so the
                 partitions depend on each other), xmin and xmax reflecting, corner source.
  nonconvex:     shestakov_quad_a0.05_n8.vtu (cycles within a partition), xmin and xmax
                 reflecting, corner source.
Each runs with classic Richardson and with GMRES.
  ADJ_DSA_MAX_REL_DIFF: max |(-1)^ell phi_adj - phi_fwd^T| / max |phi_fwd^T| over all cases (0 to
    roundoff; forward diffusion coefficients in adjoint mode gave about 1e-3).
  ADJ_DSA_SWEEP_MISMATCH: the number of cases with different sweep counts (0).
The all-reflecting case is also converged and compared with the exact infinite-medium solution
phi = (diag(sigma_t) - S0^T)^-1 q.
  ADJ_DSA_INF_MEDIUM_REL_ERR: the largest relative error of the group averages (0 to the solver
    tolerance 1e-12).
With use_gpus=True the sweeps run on the GPU.
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
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator, LinearGraphPartitioner, OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier

    def global_max(value):
        return MPIAllReduce(value, "max")


SIGMA_T = (1.0, 2.0)
# S[ell][to][from], nonsymmetric, with up- and downscatter in both moments.
S = (((0.60, 0.05), (0.38, 1.90)),
     ((0.20, 0.00), (0.15, 0.60)))
Q = (1.0, 0.5)
MAX_ITS = 5
NUM_GROUPS = 2
NUM_MOMENTS = 3


def get_option(name, default):
    return globals()[name] if name in globals() else default


def write_xs(file_name, transpose):
    if rank == 0:
        lines = ["NUM_GROUPS 2", "NUM_MOMENTS 2", "SIGMA_T_BEGIN"]
        lines += [f"{g} {SIGMA_T[g]}" for g in range(2)]
        lines += ["SIGMA_T_END", "TRANSFER_MOMENTS_BEGIN"]
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
    barrier()
    if rank == 0:
        os.remove(file_name)
    return xs


def ortho(nx, ny, linear=False):
    xs_ = [i / nx for i in range(nx + 1)]
    ys_ = [i / ny for i in range(ny + 1)]
    kwargs = {"partitioner": LinearGraphPartitioner()} if linear else {}
    return OrthogonalMeshGenerator(node_sets=[xs_, ys_], **kwargs).Execute()


def nonconvex():
    grid = FromFileMeshGenerator(
        filename="../../../../assets/mesh/shestakov_quad_a0.05_n8.vtu"
    ).Execute()
    grid.SetOrthogonalBoundaries()
    return grid


def make_problem(grid, xs, adjoint, sources, reflecting, method, max_its, tol):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=4, n_azimuthal=8, scattering_order=1
                ),
                "inner_linear_method": method,
                "l_abs_tol": tol,
                "l_max_its": max_its,
                "apply_wgdsa": True,
                "wgdsa_l_abs_tol": 1.0e-13,
                "apply_tgdsa": True,
                "tgdsa_l_abs_tol": 1.0e-13,
            }
        ],
        xs_map=[{"block_ids": [0, 1], "xs": xs}],
        volumetric_sources=sources,
        boundary_conditions=[{"name": n, "type": "reflecting"} for n in reflecting],
        options={"adjoint": adjoint, "verbose_inner_iterations": False},
        **({"use_gpus": True} if get_option("use_gpus", False) else {}),
    )


def moment_sign(k):
    # Local flux moments are stored by node, then moment, then group. The 2D P1 moments are (0, 0),
    # (1, -1), and (1, 1). The physical adjoint is the forward solution with S^T reflected, so its
    # moments of degree ell carry the factor (-1)^ell.
    return 1.0 if (k // NUM_GROUPS) % NUM_MOMENTS == 0 else -1.0


def solve(problem):
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return solver.GetNumSweeps(), list(problem.GetPhiNewLocal())


def group_avg(problem, g):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("avg")
    ffi.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    ffi.AddFieldFunction(problem.GetScalarFluxFieldFunction()[g])
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    xs_s = write_xs("adjoint_dsa_equivalence_S.xs", transpose=False)
    xs_st = write_xs("adjoint_dsa_equivalence_ST.xs", transpose=True)
    corner = RPPLogicalVolume(xmin=0.0, xmax=0.3, ymin=0.0, ymax=0.3, infz=True)

    def uniform_source():
        return [VolumetricSource(block_ids=[0, 1], group_strength=list(Q))]

    def corner_source():
        return [VolumetricSource(block_ids=[1], group_strength=list(Q))]

    all_faces = ("xmin", "xmax", "ymin", "ymax")
    cases = (
        ("inf_xy", lambda: ortho(10, 10), uniform_source, all_faces),
        ("cut_row_mixed", lambda: ortho(10, 9, linear=True), corner_source, ("xmin", "xmax")),
        ("nonconvex", nonconvex, corner_source, ("xmin", "xmax")),
    )

    max_rel_diff = 0.0
    sweep_mismatch = 0
    for name, make_grid, make_source, reflecting in cases:
        grid = make_grid()
        grid.SetUniformBlockID(0)
        grid.SetBlockIDFromLogicalVolume(corner, 1, True)
        for method in ("classic_richardson", "petsc_gmres"):
            sweeps_adj, phi_adj = solve(
                make_problem(grid, xs_s, True, make_source(), reflecting, method, MAX_ITS, 1e-14)
            )
            sweeps_fwd, phi_fwd = solve(
                make_problem(grid, xs_st, False, make_source(), reflecting, method, MAX_ITS, 1e-14)
            )
            scale = global_max(max(abs(v) for v in phi_fwd))
            diff = global_max(
                max(abs(moment_sign(k) * a - b) for k, (a, b) in enumerate(zip(phi_adj, phi_fwd)))
            )
            max_rel_diff = max(max_rel_diff, diff / scale)
            sweep_mismatch += int(sweeps_adj != sweeps_fwd)

    # Converged infinite medium: phi = (diag(sigma_t) - S0^T)^-1 q.
    a00, a01 = SIGMA_T[0] - S[0][0][0], -S[0][1][0]
    a10, a11 = -S[0][0][1], SIGMA_T[1] - S[0][1][1]
    det = a00 * a11 - a01 * a10
    exact = ((a11 * Q[0] - a01 * Q[1]) / det, (a00 * Q[1] - a10 * Q[0]) / det)
    grid = ortho(10, 10)
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(corner, 1, True)
    inf_err = 0.0
    for method in ("classic_richardson", "petsc_gmres"):
        problem = make_problem(grid, xs_s, True, uniform_source(), all_faces, method, 1000, 1e-12)
        solve(problem)
        for g in range(2):
            inf_err = max(inf_err, abs(group_avg(problem, g) / exact[g] - 1.0))

    if rank == 0:
        print(f"ADJ_DSA_MAX_REL_DIFF={max_rel_diff:.6e}")
        print(f"ADJ_DSA_SWEEP_MISMATCH={sweep_mismatch}")
        print(f"ADJ_DSA_INF_MEDIUM_REL_ERR={inf_err:.6e}")
