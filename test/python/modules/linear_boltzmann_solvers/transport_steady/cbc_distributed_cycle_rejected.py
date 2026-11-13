#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT
"""A two-cell cross-rank cycle must fail during CBC setup on both ranks.

Two faces connect the same two cells with opposite upwind relationships for
an oblique ordinate. Each rank owns one cell, hence neither local graph cycles.
The hand constructed geometry has no numerical reference: this is an exact
structural rejection test, registered with a bounded runner timeout.
"""
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank, size = MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.xs import MultiGroupXS
    from pyopensn.solver import DiscreteOrdinatesProblem

if __name__ == "__main__":
    if size != 2:
        sys.exit("This test requires two ranks")
    grid = ExtruderMeshGenerator(
        inputs=[FromFileMeshGenerator(filename="cbc_distributed_cycle.obj")],
        layers=[{"z": 0.25, "n": 1}],
        partitioner=KBAGraphPartitioner(ny=2, ycuts=[0.6]),
    ).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)
    try:
        DiscreteOrdinatesProblem(
            mesh=grid, num_groups=1,
            groupsets=[{"groups_from_to": (0, 0), "allow_cycles": False,
                        "angular_quadrature": GLCProductQuadrature3DXYZ(
                            n_polar=2, n_azimuthal=8, scattering_order=0)}],
            xs_map=[{"block_ids": [0], "xs": xs}], sweep_type="CBC",
        )
    except RuntimeError as error:
        message = str(error)
        if not all(token in message for token in
                   ("CBC", "ordering", "global cell dependency cycle", "rank", "cell")):
            raise
        print(f"CBC_GLOBAL_CYCLE_REJECTED_RANK_{rank}=1", flush=True)
        print(message, flush=True)
    else:
        sys.exit("CBC accepted a genuine cross-rank cell cycle")
