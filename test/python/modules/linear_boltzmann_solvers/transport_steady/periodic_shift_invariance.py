#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""A half-period source shift must produce the same shift in the slab flux."""

import os
import sys

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def solve(source_start, sweep_type):
    n_cells = 20
    nodes = [i / n_cells for i in range(n_cells + 1)]
    mesh = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    mesh.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)
    source_region = RPPLogicalVolume(
        infx=True, infy=True, zmin=source_start, zmax=source_start + 0.5
    )
    source = VolumetricSource(logical_volume=source_region, group_strength=[1.0])
    quadrature = GLProductQuadrature1DSlab(n_polar=8, scattering_order=0)
    problem = DiscreteOrdinatesProblem(
        mesh=mesh,
        num_groups=1,
        sweep_type=sweep_type,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-11,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source],
        boundary_conditions=[
            {"name": "zmin", "type": "periodic"},
            {"name": "zmax", "type": "periodic"},
        ],
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    return list(problem.GetPhiNewLocal()), solver.ComputeBalanceTable()


if __name__ == "__main__":
    if size != 1:
        sys.exit("Shift-invariance test requires one MPI rank.")
    for sweep in ("AAH", "CBC"):
        base, base_balance = solve(0.0, sweep)
        shifted, shifted_balance = solve(0.5, sweep)
        half = len(base) // 2
        error = max(abs(base[i] - shifted[(i + half) % len(base)]) for i in range(len(base)))
        contrast = max(base) - min(base)
        print(f"PeriodicShiftError{sweep}={error:.8e}")
        print(f"PeriodicShiftContrast{sweep}={contrast:.8e}")
        print(f"PeriodicShiftBalance{sweep}={base_balance['balance']:.8e}")
        print(f"PeriodicShiftedBalance{sweep}={shifted_balance['balance']:.8e}")
