#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT
"""A local cell cycle on rank 0 must reject CBC setup on both ranks.

The cyclic3 mesh with an x cut at 0.8 leaves only an acyclic outer strip
on rank 1. This checks collective rejection, not a transport solution.
"""
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.xs import MultiGroupXS
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

if __name__ == "__main__":
    if size != 2:
        sys.exit("This test requires two ranks.")
    grid = ExtruderMeshGenerator(
        inputs=[FromFileMeshGenerator(
            filename="../../../../assets/mesh/square2x2_partition_cyclic3.obj")],
        layers=[{"z": 1.0, "n": 1}],
        partitioner=KBAGraphPartitioner(nx=2, xcuts=[0.8]),
    ).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)
    quadrature = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=8, scattering_order=0)
    try:
        problem = DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[{
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "allow_cycles": False,
            }],
            xs_map=[{"block_ids": [0], "xs": xs}],
            sweep_type="CBC",
        )
        SteadyStateSourceSolver(problem=problem).Initialize()
    except RuntimeError as error:
        expected = ("CBC sweep setup failed on rank 0: CBC_SPDS: "
                    "Cyclic dependencies found in the local cell graph.")
        if not str(error).startswith(expected):
            raise
        print(f"CBC_LOCAL_CYCLE_REJECTED_RANK_{rank}=1", flush=True)
        print(f"Rank {rank}: {error}", flush=True)
    else:
        sys.exit("CBC accepted a local cell cycle with allow_cycles=False.")
