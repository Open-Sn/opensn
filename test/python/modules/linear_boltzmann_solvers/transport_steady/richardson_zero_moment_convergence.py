#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""Classic Richardson must converge when higher flux moments are zero.

A homogeneous slab with reflecting boundaries and a uniform isotropic source is an infinite
medium: the angular flux is isotropic and uniform, so every flux moment above the scalar flux is
zero and holds only round-off. The pointwise convergence test must judge those moments against the
local scalar flux. Judged against their own round-off magnitude, they keep changing until the
iterate is stationary to machine precision, which takes about ln(eps) / ln(c) ~ 340 iterations,
whatever the tolerance. With tolerance 1e-6 and c = 0.9, the iteration needs about
ln(1e-6 (1 - c)) / ln(c) ~ 150 iterations, so the solve must converge within its 250-iteration
limit ("status = converged").

The scalar flux must equal the infinite-medium solution q / (sigma_t (1 - c)) within twice the
tolerance: the stopping test estimates the error as change / (1 - rho), which is not a bound (the
relative error at convergence is 8.6e-7).
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":
    sigma_t = 1.0
    c = 0.9
    q = 1.0

    grid = OrthogonalMeshGenerator(node_sets=[[0.25 * i for i in range(9)]]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=c)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=7),
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 250,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[q])],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    phi = problem.GetScalarFluxFieldFunction()[0]
    exact = q / (sigma_t * (1.0 - c))
    errors = []
    for op in ("min", "max"):
        ffi = FieldFunctionInterpolationVolume()
        ffi.SetOperationType(op)
        ffi.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
        ffi.AddFieldFunction(phi)
        ffi.Execute()
        errors.append(abs(ffi.GetValue() - exact) / exact)
    if rank == 0:
        print(f"PHI_MAX_REL_ERR={max(errors):.6e}")
