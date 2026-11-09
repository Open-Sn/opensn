#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Source moments made from an adjoint flux are in the solver orientation that source moments use.

An adjoint solve leaves its flux moments in the physical orientation, in which a moment of degree
ell carries a factor (-1)^ell. Source moments are applied in the solver orientation, so moments made
from a physical adjoint flux must be mapped back. The check is an exact linear identity. With L0 the
streaming and collision operator and S^T the transposed scattering operator, the adjoint flux of
problem A satisfies

    phi_A = L0^-1 q + L0^-1 (S^T phi_A).

  A  adjoint, two groups with up- and downscatter, P1 scattering, volumetric source q. After the
     solve, CreateAndWriteSourceMoments writes s = S^T phi_A and WriteFluxMoments writes phi_A.
  C  adjoint, same sigma_t, no scattering, use_source_moments with s read back: L0^-1 s.
  E  adjoint with A's cross sections and use_source_moments. ReadFluxMomentsAndMakeSourceMoments
     makes s from the phi_A file, and the cross sections are then replaced by the absorber, so the
     solve gives L0^-1 s as well.
  D  adjoint, same sigma_t, no scattering, volumetric source q: L0^-1 q.

So phi_C + phi_D = phi_A and phi_E + phi_D = phi_A in every moment, to the solver tolerance (1e-13).
With wrongly oriented moments the P1 part of s has the wrong sign, which gives errors of a few
percent.
  ADJ_SRC_MOMENTS_CREATE_ERR, ADJ_SRC_MOMENTS_READ_ERR: max |phi_A - phi_X - phi_D| / max |phi_A|
    for X = C and X = E.
  ADJ_SRC_MOMENTS_FORWARD_FILE_REJECTED: 1 when a forward problem rejects the adjoint flux file in
    ReadFluxMomentsAndMakeSourceMoments.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank

    def global_max(value):
        return MPI.COMM_WORLD.allreduce(value, op=MPI.MAX)

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:

    def global_max(value):
        return MPIAllReduce(value, "max")


SIGMA_T = [1.0, 1.2]
# (moment, from group, to group, value)
TRANSFERS = [
    (0, 0, 0, 0.3), (0, 0, 1, 0.3), (0, 1, 0, 0.05), (0, 1, 1, 0.5),
    (1, 0, 0, 0.15), (1, 0, 1, 0.1), (1, 1, 0, 0.02), (1, 1, 1, 0.2),
]


def make_xs(transfers, tag):
    lines = ["NUM_GROUPS 2", "NUM_MOMENTS 2", "SIGMA_T_BEGIN"]
    lines += [f"{g} {SIGMA_T[g]}" for g in range(2)]
    lines += ["SIGMA_T_END", "TRANSFER_MOMENTS_BEGIN"]
    lines += [f"M_GFROM_GTO_VAL {m} {gf} {gt} {v}" for m, gf, gt, v in transfers]
    lines += ["TRANSFER_MOMENTS_END"]
    path = f"adjoint_source_moments_{tag}_{rank}.xs"
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    os.remove(path)
    return xs


if __name__ == "__main__":
    nodes = [i * 0.1 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(infx=True, infy=True, zmin=1.0, zmax=2.01), 1, True
    )
    scatterer = make_xs(TRANSFERS, "scatterer")
    absorber = make_xs([(0, 0, 0, 0.0)], "absorber")

    def make_problem(xs, adjoint=True, use_moments=False, source=True):
        sources = [VolumetricSource(block_ids=[1], group_strength=[0.0, 1.0])] if source else []
        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=2,
            groupsets=[
                {
                    "groups_from_to": (0, 1),
                    "angular_quadrature": GLProductQuadrature1DSlab(
                        n_polar=16, scattering_order=1
                    ),
                    "inner_linear_method": "petsc_gmres",
                    "l_abs_tol": 1.0e-13,
                    "l_max_its": 300,
                }
            ],
            xs_map=[{"block_ids": [0, 1], "xs": xs}],
            volumetric_sources=sources,
            options={
                "adjoint": adjoint,
                "use_source_moments": use_moments,
                "verbose_inner_iterations": False,
            },
        )

    def solve(problem):
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()
        return list(problem.GetPhiNewLocal())

    problem_a = make_problem(scatterer)
    phi_a = solve(problem_a)
    problem_a.CreateAndWriteSourceMoments("adjoint_source_moments_s")
    problem_a.WriteFluxMoments("adjoint_source_moments_phi")

    problem_c = make_problem(absorber, use_moments=True, source=False)
    problem_c.ReadSourceMoments("adjoint_source_moments_s", False)
    phi_c = solve(problem_c)

    problem_e = make_problem(scatterer, use_moments=True, source=False)
    problem_e.ReadFluxMomentsAndMakeSourceMoments("adjoint_source_moments_phi", False)
    problem_e.SetXSMap(xs_map=[{"block_ids": [0, 1], "xs": absorber}])
    phi_e = solve(problem_e)

    phi_d = solve(make_problem(absorber))

    forward = make_problem(scatterer, adjoint=False, use_moments=True, source=False)
    try:
        forward.ReadFluxMomentsAndMakeSourceMoments("adjoint_source_moments_phi", False)
        forward_file_rejected = 0
    except RuntimeError as error:
        forward_file_rejected = int("does not match the requested role" in str(error))

    def rel_err(phi_x):
        local = max(abs(a - x - d) for a, x, d in zip(phi_a, phi_x, phi_d))
        scale = max(abs(a) for a in phi_a)
        return global_max(local) / global_max(scale)

    create_err = rel_err(phi_c)
    read_err = rel_err(phi_e)
    for prefix in ("adjoint_source_moments_s", "adjoint_source_moments_phi"):
        os.remove(f"{prefix}{rank}.h5")
    if rank == 0:
        print(f"ADJ_SRC_MOMENTS_CREATE_ERR={create_err:.6e}")
        print(f"ADJ_SRC_MOMENTS_READ_ERR={read_err:.6e}")
        print(f"ADJ_SRC_MOMENTS_FORWARD_FILE_REJECTED={forward_file_rejected}")
