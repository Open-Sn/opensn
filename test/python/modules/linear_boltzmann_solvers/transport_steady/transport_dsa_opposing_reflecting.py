#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
WGDSA with opposing reflecting boundaries.

An infinite medium, represented by a 2D square and a 3D cube with reflecting
boundaries on every face, holds a 1-group material with sigma_t = 1 cm^-1 and
scattering ratio c = 0.99 and a uniform source Q = 1. The exact solution is
phi = Q / (sigma_t (1 - c)) = 100.

The incoming angular flux on opposing reflecting boundaries is lagged by one
iteration and carried as an unknown of the groupset iteration. The WGDSA
correction must be applied to these lagged fluxes as well; without it, source
iteration (Richardson) with WGDSA diverges on this problem. Each inner solver
is given 60 iterations, which unaccelerated source iteration needs thousands of
to converge (spectral radius ~0.99).

Check:
  DSA_REFLECTING_MAX_REL_ERR: max over dimensions and inner solvers
  (classic_richardson, petsc_richardson, petsc_gmres) of |phi / 100 - 1|.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY, GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def field_stat(problem, op):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType(op)
    ffi.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    ffi.AddFieldFunction(problem.GetScalarFluxFieldFunction()[0])
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    sigma_t, c, Q = 1.0, 0.99, 1.0
    phi_exact = Q / (sigma_t * (1.0 - c))

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, c)

    max_err = 0.0
    for dim in (2, 3):
        n = 8 if dim == 2 else 4
        nodes = [i * 10.0 / n for i in range(n + 1)]
        grid = OrthogonalMeshGenerator(node_sets=[nodes] * dim).Execute()
        grid.SetUniformBlockID(0)
        if dim == 2:
            quad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8, scattering_order=0)
            names = ("xmin", "xmax", "ymin", "ymax")
        else:
            quad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)
            names = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")

        for method in ("classic_richardson", "petsc_richardson", "petsc_gmres"):
            problem = DiscreteOrdinatesProblem(
                mesh=grid,
                num_groups=1,
                groupsets=[
                    {
                        "groups_from_to": (0, 0),
                        "angular_quadrature": quad,
                        "inner_linear_method": method,
                        "l_abs_tol": 1.0e-10,
                        "l_max_its": 60,
                        "apply_wgdsa": True,
                        "wgdsa_l_abs_tol": 1.0e-12,
                    }
                ],
                xs_map=[{"block_ids": [0], "xs": xs}],
                volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[Q])],
                boundary_conditions=[{"name": name, "type": "reflecting"} for name in names],
                options={"verbose_inner_iterations": False},
            )
            solver = SteadyStateSourceSolver(problem=problem)
            solver.Initialize()
            solver.Execute()
            for op in ("min", "max"):
                err = abs(field_stat(problem, op) / phi_exact - 1.0)
                max_err = max(max_err, err if err == err else 1.0e300)

    if rank == 0:
        print(f"DSA_REFLECTING_MAX_REL_ERR {max_err:.6e}")
