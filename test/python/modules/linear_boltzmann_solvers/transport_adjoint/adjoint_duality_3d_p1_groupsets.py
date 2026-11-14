#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Forward-adjoint duality for a 3D, two-group problem with P1 scattering in two groupsets (2 ranks).

The discrete transport operator L satisfies <q_adj, L^-1 q> = <(L^T)^-1 q_adj, q>, so the forward
response of a detector equals the adjoint flux folded with the forward source:

    R_fwd = sum_g int_B sigma_det,g phi_g dV = sum_g int_A q_g phi_adj,g dV = R_adj.

The unit cube has 5x5x5 cells. A forward source q = (1, 0.5) fills block 1, [0, 0.4]^3, and the
detector sigma_det = (0.3, 0.8) fills block 2, [0.6, 1]^3. xmin is reflecting and the other faces
are vacuum. The scattering has up- and downscatter in both the P0 and the P1 moment, and each group
is its own groupset, so the identity tests the transposition of every moment across groupsets. The
forward and adjoint problems are solved independently to 1e-12, and both integrals are evaluated by
the response evaluator with the consistent mass matrix.
  ADJ_DUALITY_3D_REL_DIFF: |R_adj - R_fwd| / |R_fwd| (0 to the solver tolerance).
  ADJ_DUALITY_3D_RESPONSE: R_fwd, which shows the detector sees the source.
The sweep_type option selects AAH (default) or CBC sweeps, and use_gpus=True runs them on the GPU.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.response import ResponseEvaluator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier

SIGMA_T = (1.0, 1.5)
# (moment, from group, to group, value)
TRANSFERS = (
    (0, 0, 0, 0.40), (0, 0, 1, 0.30), (0, 1, 0, 0.08), (0, 1, 1, 0.90),
    (1, 0, 0, 0.15), (1, 0, 1, 0.10), (1, 1, 0, 0.03), (1, 1, 1, 0.30),
)
Q = [1.0, 0.5]
SIGMA_DET = [0.3, 0.8]
SWEEP_TYPE = globals().get("sweep_type", "AAH")
USE_GPUS = globals().get("use_gpus", False)


def load_xs(file_name):
    if rank == 0:
        lines = ["NUM_GROUPS 2", "NUM_MOMENTS 2", "SIGMA_T_BEGIN"]
        lines += [f"{g} {SIGMA_T[g]}" for g in range(2)]
        lines += ["SIGMA_T_END", "TRANSFER_MOMENTS_BEGIN"]
        lines += [f"M_GFROM_GTO_VAL {m} {gf} {gt} {v}" for m, gf, gt, v in TRANSFERS]
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


if __name__ == "__main__":
    nodes = [i * 0.2 for i in range(6)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.0, xmax=0.4, ymin=0.0, ymax=0.4, zmin=0.0, zmax=0.4), 1, True
    )
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.6, xmax=1.0, ymin=0.6, ymax=1.0, zmin=0.6, zmax=1.0), 2, True
    )
    xs = load_xs("adjoint_duality_3d_p1_groupsets.xs")

    def make_problem(adjoint, block, strength):
        quadrature = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=1)
        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=2,
            groupsets=[
                {
                    "groups_from_to": (g, g),
                    "angular_quadrature": quadrature,
                    "inner_linear_method": "petsc_gmres",
                    "l_abs_tol": 1.0e-13,
                    "l_max_its": 300,
                }
                for g in range(2)
            ],
            xs_map=[{"block_ids": [0, 1, 2], "xs": xs}],
            volumetric_sources=[VolumetricSource(block_ids=[block], group_strength=strength)],
            boundary_conditions=[{"name": "xmin", "type": "reflecting"}],
            sweep_type=SWEEP_TYPE,
            options={
                "adjoint": adjoint,
                "max_ags_iterations": 500,
                "ags_tolerance": 1.0e-12,
                "verbose_inner_iterations": False,
                "verbose_outer_iterations": False,
            },
            **({"use_gpus": True} if USE_GPUS else {}),
        )

    def response(adjoint, source_block, source_strength, detector_block, detector_strength):
        problem = make_problem(adjoint, source_block, source_strength)
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()
        prefix = f"adjoint_duality_3d_p1_groupsets_{'adj' if adjoint else 'fwd'}"
        problem.WriteFluxMoments(prefix)
        evaluator = ResponseEvaluator(problem=problem)
        evaluator.SetOptions(
            data=[{"name": "phi", "file_prefixes": {"flux_moments": prefix}}],
            sources={"material": [{"block_id": detector_block, "strength": detector_strength}]},
        )
        value = evaluator.EvaluateResponse("phi")
        os.remove(f"{prefix}{rank}.h5")
        return value

    r_fwd = response(False, 1, Q, 2, SIGMA_DET)
    r_adj = response(True, 2, SIGMA_DET, 1, Q)

    if rank == 0:
        print(f"ADJ_DUALITY_3D_REL_DIFF={abs(r_adj - r_fwd) / abs(r_fwd):.6e}")
        print(f"ADJ_DUALITY_3D_RESPONSE={r_fwd:.10e}")
