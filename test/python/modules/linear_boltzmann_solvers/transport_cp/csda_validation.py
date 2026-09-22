#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""CSDA restrictions that must hold after construction as well as during it.

Checks that a CSDA problem rejects an uncollided flux file at construction,
rejects a later switch to adjoint mode, and rejects replacement cross sections
whose charged-particle block would be split across groupsets or whose group
count does not match the problem. Supplied material energy structures must agree,
and at least one is required. All rejections must be ValueErrors. Each rejected
call must leave the problem unchanged, so the original problem must still
solve to a converged particle and energy balance afterward.
"""

import os
import struct
import tempfile

from mpi4py import MPI

if "opensn_console" not in globals():
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

comm = MPI.COMM_WORLD


def make_xs(stopping, bounds=None):
    """Write a CEPXS CSDA fixture with len(stopping) groups, shared by ranks."""
    num_groups = len(stopping)
    path = None
    if comm.rank == 0:
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), suffix=".bxslib", delete=False) as f:
            path = f.name

            def record(payload):
                marker = struct.pack("<I", len(payload))
                f.write(marker + payload + marker)

            table = [0.0] * (11 * num_groups)
            for g in range(num_groups):
                table[11 * g + 1] = [0.05, 0.03, 0.01, 0.005][g]
                table[11 * g + 2] = [0.5, 0.35, 0.2, 0.1][g]
                table[11 * g + 4] = stopping[g]
                table[11 * g + 7] = [1.0, 1.2, 1.5, 1.8][g]
            if bounds is None:
                bounds = [float(num_groups - i) for i in range(num_groups + 1)]
            record(b"CSDA VALIDATION REGRESSION")
            record(struct.pack("<8i", num_groups, 1, 11, 8, 9, 1, 0, 1))
            record(struct.pack(f"<{num_groups + 1}d", *bounds))
            record(struct.pack(f"<{11 * num_groups}d", *table))
    path = comm.bcast(path, root=0)
    xs = MultiGroupXS()
    xs.LoadFromCEPXS(path, material_id=0, csda_format=True)
    comm.Barrier()
    if comm.rank == 0:
        os.remove(path)
    return xs


def expect_rejected(action, message, label):
    """Require a ValueError, as documented, whose text contains message."""
    try:
        action()
    except ValueError as err:
        assert message in str(err), (label, str(err))
        return
    raise AssertionError(f"{label}: expected ValueError containing '{message}'")


grid = OrthogonalMeshGenerator(node_sets=[[i / 8 for i in range(9)]]).Execute()
grid.SetUniformBlockID(0)

# Groups 1 and 2 are charged, and both sit in the second groupset.
split_safe_xs = make_xs([0.0, 0.3, 0.2])
# Groups 0-2 are all charged, so this block would span both groupsets.
split_xs = make_xs([0.4, 0.3, 0.2])
# A fourth group: more groups than the problem, which CSDA does not allow.
four_group_xs = make_xs([0.0, 0.3, 0.2, 0.1])


def problem_args(**overrides):
    quad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)
    args = dict(
        mesh=grid,
        num_groups=3,
        groupsets=[
            {"groups_from_to": (0, 0), "angular_quadrature": quad,
             "inner_linear_method": "petsc_gmres", "l_abs_tol": 1e-12, "l_max_its": 1000},
            {"groups_from_to": (1, 2), "angular_quadrature": quad,
             "inner_linear_method": "petsc_gmres", "l_abs_tol": 1e-12, "l_max_its": 1000},
        ],
        xs_map=[{"block_ids": [0], "xs": split_safe_xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[0.0, 2.0, 0.0])],
        boundary_conditions=[{"name": name, "type": "reflecting"} for name in ("zmin", "zmax")],
        options={"csda_enabled": True},
    )
    args.update(overrides)
    return args


# An uncollided flux file is rejected at construction, before the file is read.
expect_rejected(
    lambda: DiscreteOrdinatesProblem(**problem_args(uncollided_flux="does_not_exist.h5")),
    "CSDA is not supported with an uncollided flux file",
    "uncollided flux")

# CSDA requires at least one material to supply energy bounds.
unbounded = MultiGroupXS()
unbounded.CreateSimpleOneGroup(sigma_t=0.0, c=0.0)
expect_rejected(
    lambda: DiscreteOrdinatesProblem(**problem_args(
        num_groups=1,
        groupsets=[{"groups_from_to": (0, 0),
                    "angular_quadrature": GLProductQuadrature1DSlab(
                        n_polar=4, scattering_order=0)}],
        xs_map=[{"block_ids": [0], "xs": unbounded}], volumetric_sources=[])),
    "complete energy group structure", "missing problem energy structure")

# CSDA with no charged group would silently run neutral transport, so it is rejected.
neutral_xs = make_xs([0.0, 0.0, 0.0])
expect_rejected(
    lambda: DiscreteOrdinatesProblem(**problem_args(
        xs_map=[{"block_ids": [0], "xs": neutral_xs}], volumetric_sources=[])),
    "at least one charged-particle group", "no charged groups")

problem = DiscreteOrdinatesProblem(**problem_args())

expect_rejected(lambda: problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": neutral_xs}]),
                "at least one charged-particle group", "replacement with no charged groups")

# Switching an existing CSDA problem to adjoint mode is rejected.
expect_rejected(lambda: problem.SetAdjoint(True),
                "CSDA is not supported for adjoint problems", "adjoint switch")

# Replacement cross sections that split a charged block across groupsets are rejected.
expect_rejected(lambda: problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": split_xs}]),
                "is split across groupsets", "split XS map")

# Replacement CSDA cross sections with the wrong number of groups are rejected as bad
# input rather than as an internal error.
expect_rejected(lambda: problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": four_group_xs}]),
                "incompatible with the configured number of groups", "group count")

# Even a neutral material must agree on both bounds, not just group widths.
for bounds in ([6.0, 4.0, 2.0, 0.0], [4.0, 3.0, 2.0, 1.0]):
    conflicting_xs = make_xs([0.0, 0.0, 0.0], bounds=bounds)
    conflicting_map = [{"block_ids": [0], "xs": split_safe_xs},
                       {"block_ids": [1], "xs": conflicting_xs}]
    expect_rejected(
        lambda: DiscreteOrdinatesProblem(**problem_args(xs_map=conflicting_map)),
        "all supplied energy group structures must match", "constructor energy mismatch")
    expect_rejected(
        lambda: problem.SetXSMap(xs_map=conflicting_map),
        "all supplied energy group structures must match",
        "replacement energy mismatch")

# The rejected calls must not have changed the problem: it still solves in forward
# mode with the original cross sections and closes both balances.
solver = SteadyStateSourceSolver(problem=problem)
solver.Initialize()
solver.Execute()
balance = solver.ComputeBalanceTable()
max_balance = max(abs(balance["csda_particle_balance"]), abs(balance["csda_energy_balance"]))
assert max_balance < 1e-9, balance

if comm.rank == 0:
    print(f"CSDA_VALIDATION_MAX_BALANCE={max_balance:.12e}")
