#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""One-group periodic transport in Cartesian slab, XY, and XYZ meshes.

With a spatially uniform isotropic source Q=1 and pure absorption Sigma_a=1,
the exact infinite-medium solution is phi=1 everywhere. Fully periodic
boundaries must reproduce it independently of mesh dimension or partition.
"""

import os
import sys

from mpi4py import MPI

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import (
        GLProductQuadrature1DSlab,
        GLCProductQuadrature2DXY,
        GLCProductQuadrature3DXYZ,
    )
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def solve(dimension, sweep_type):
    nodes = [0.0, 0.25, 0.5, 0.75, 1.0] if dimension == 1 else [0.0, 0.5, 1.0]
    node_sets = [nodes] * dimension
    mesh = OrthogonalMeshGenerator(node_sets=node_sets).Execute()
    mesh.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)
    source = VolumetricSource(block_ids=[0], group_strength=[1.0])

    if dimension == 1:
        quadrature = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)
        names = ("zmin", "zmax")
    elif dimension == 2:
        quadrature = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)
        names = ("xmin", "xmax", "ymin", "ymax")
    else:
        quadrature = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)
        names = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")

    problem = DiscreteOrdinatesProblem(
        mesh=mesh,
        num_groups=1,
        sweep_type=sweep_type,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source],
        boundary_conditions=[{"name": name, "type": "periodic"} for name in names],
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    local_phi = list(problem.GetPhiNewLocal())
    local_error = max((abs(phi - 1.0) for phi in local_phi), default=0.0)
    error = MPI.COMM_WORLD.allreduce(local_error, op=MPI.MAX)
    balance = solver.ComputeBalanceTable()
    if MPI.COMM_WORLD.rank == 0:
        print(f"PeriodicError{dimension}D{sweep_type}={error:.8e}")
        print(f"PeriodicBalance{dimension}D{sweep_type}={balance['balance']:.8e}")


if __name__ == "__main__":
    for dim in (1, 2, 3):
        for sweep in ("AAH", "CBC"):
            solve(dim, sweep)
