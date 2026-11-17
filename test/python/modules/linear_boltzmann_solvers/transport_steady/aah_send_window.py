#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""
Limiting the posted AAH sweep sends with ``max_mpi_outstanding_sends`` must not change the solution.

The mesh is the extruded ``square2x2_partition_cyclic3`` mesh on 4 ranks. Its partitions form
cycles in the sweep graph between ranks, so some messages go to delayed successors, which take them
only after the sweep. The limited run allows 1 posted send per rank and splits the messages at
``max_mpi_message_size = 1024`` bytes, so most angle sets have more messages than the limit and post
them one at a time. With CSDA a psi message and its psiE partner form one unit of 2 sends, which is
larger than the limit and is admitted only when no send is outstanding.

Claim: the limited and the unlimited runs agree bitwise. Both use ``deterministic_sweeps``, so the
moments and outflows do not depend on the order in which messages arrive, and the messages carry the
same doubles. Agreement alone does not show that either is correct, so the limited run must also
conserve particles: production equals absorption plus outflow, relative to production (with the
CSDA deposition out of the lowest group in the CSDA case). GMRES is converged to an absolute
tolerance of 1e-12, so the balance closes to about that level. The option must also be rejected for
CBC sweeps, which do not limit their sends. Sending the split messages from the wrong buffer offset
breaks both the bitwise agreement and the balance.

Cases: ``steady`` (10-group fixed source that uses every group of its cross sections, so the energy
range is closed) and ``csda`` (3-group CSDA fixed source).
"""

import hashlib
import os
import struct
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    comm = MPI.COMM_WORLD

    def MPIAllReduce(value, op="sum"):
        return comm.allreduce(value, op={"sum": MPI.SUM, "min": MPI.MIN}[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


CASE = globals().get("case", "steady")


def make_grid():
    grid = ExtruderMeshGenerator(
        inputs=[
            FromFileMeshGenerator(
                filename="../../../../assets/mesh/square2x2_partition_cyclic3.obj"
            )
        ],
        layers=[{"z": 0.4, "n": 2}, {"z": 0.8, "n": 2}],
        partitioner=KBAGraphPartitioner(nx=2, ny=2, xcuts=[0.0], ycuts=[0.0]),
    ).Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    return grid


def load_xs():
    xs = MultiGroupXS()
    if CASE == "csda":
        xs.LoadFromCEPXS(
            "../../../../assets/xs/cepxs_synthetic_csda_3g.bxslib", material_id=0, csda_format=True
        )
    else:
        xs.LoadFromOpenSn("../../../../assets/xs/response_2d_3_mat1.xs")
    return xs


def make_problem(grid, xs, sweep_type="AAH", **options):
    num_groups = 3 if CASE == "csda" else 10
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=4, n_azimuthal=8, scattering_order=0
                ),
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 300,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[
            VolumetricSource(block_ids=[0], group_strength=[1.0] + [0.0] * (num_groups - 1))
        ],
        options={
            "deterministic_sweeps": True,
            "csda_enabled": CASE == "csda",
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            **options,
        },
        sweep_type=sweep_type,
    )


def solve(grid, xs, **options):
    problem = make_problem(grid, xs, **options)
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    table = solver.ComputeBalanceTable()
    if CASE == "csda":
        relative_balance = table["csda_particle_balance"]
    else:
        relative_balance = table["balance"] / table["production_rate"]
    boundary_leakage = problem.ComputeLeakage([])
    leakage = [value for name in sorted(boundary_leakage) for value in boundary_leakage[name]]
    return list(problem.GetPhiNewLocal()), leakage, relative_balance


def digest(values):
    return hashlib.sha256(b"".join(struct.pack("<d", v) for v in values)).hexdigest()


def is_rejected(action):
    try:
        action()
    except ValueError as error:
        return int("max_mpi_outstanding_sends requires" in str(error))
    return 0


if __name__ == "__main__":
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")
    if CASE not in ("steady", "csda"):
        sys.exit(f"Unknown case '{CASE}'. Expected 'steady' or 'csda'.")

    grid = make_grid()
    xs = load_xs()

    unlimited = solve(grid, xs)
    limited = solve(grid, xs, max_mpi_outstanding_sends=1, max_mpi_message_size=1024)

    phi_match = MPIAllReduce(1 if digest(unlimited[0]) == digest(limited[0]) else 0, "min")
    leakage_match = 1 if digest(unlimited[1]) == digest(limited[1]) else 0
    rejected = 0
    if CASE == "steady":
        rejected = is_rejected(
            lambda: make_problem(grid, xs, sweep_type="CBC", max_mpi_outstanding_sends=8)
        )

    if rank == 0:
        print(f"SEND_WINDOW_PHI_MATCH= {phi_match}")
        print(f"SEND_WINDOW_LEAKAGE_MATCH= {leakage_match}")
        print(f"SEND_WINDOW_RELATIVE_BALANCE= {limited[2]:.6e}")
        if CASE == "steady":
            print(f"SEND_WINDOW_CBC_REJECTED= {rejected}")
