#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
WGDSA boundary residual on lagged reflecting boundaries.

An infinite medium is represented by an optically thin domain (0.1 mean free
paths across) with reflecting boundaries: a 1D slab, a 2D square (reflecting on
every face) and a 1D sphere (reflecting outer surface). It holds a 1-group
material with sigma_t = 1 cm^-1, scattering ratio c = 0.9999 and a uniform
source Q = 1. The exact solution is phi = Q / (sigma_t (1 - c)) = 1e4.

The incoming angular flux on one boundary of each pair of opposing reflecting
boundaries (and on the reflecting surface of the sphere) is lagged by one
sweep. On such a boundary the error of a sweep has an incoming part, the change
of the lagged flux, which the WGDSA correction receives as an inward current.
Without it, the correction misses the particles that the lagged boundary has not
yet returned. In this thin domain, where that boundary term dominates, the
scalar flux after four iterations is then in error by 23-59% with classic and
PETSc Richardson and by 0.06-2.4% with GMRES. With it, one accelerated
iteration removes the spatially flat error, and after four iterations the error
is below 4e-6 in every case. What remains is the slowly converging anisotropic
part of the lagged flux, which diffusion does not correct.

Each case runs exactly four iterations, so the result measures the
acceleration, not the stopping test. With classic_richardson, the boundary
residual is the change of the lagged flux over a sweep. With petsc_richardson
and petsc_gmres, WGDSA is a preconditioner applied to a residual vector whose
lagged block is the boundary residual.

Check:
  DSA_REFLECTING_RESIDUAL_MAX_REL_ERR: max over geometries and inner solvers of
  |phi / 1e4 - 1| (min and max of phi over the domain).
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
    from pyopensn.aquad import (
        GLCProductQuadrature2DXY,
        GLProductQuadrature1DSlab,
        GLProductQuadrature1DSpherical,
    )
    from pyopensn.solver import (
        DiscreteOrdinatesCurvilinearProblem,
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
    )
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
    sigma_t, c, Q, width, n = 1.0, 0.9999, 1.0, 0.1, 10
    phi_exact = Q / (sigma_t * (1.0 - c))
    nodes = [i * width / n for i in range(n + 1)]

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, c)

    cases = (
        ("slab", [nodes], GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
         ("zmin", "zmax")),
        ("xy", [nodes, nodes],
         GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8, scattering_order=0),
         ("xmin", "xmax", "ymin", "ymax")),
        ("sphere", [nodes], GLProductQuadrature1DSpherical(n_polar=8, scattering_order=0),
         ("zmax",)),
    )

    max_err = 0.0
    for geometry, node_sets, quad, names in cases:
        coord_sys = "spherical" if geometry == "sphere" else "cartesian"
        grid = OrthogonalMeshGenerator(node_sets=node_sets, coord_sys=coord_sys).Execute()
        grid.SetUniformBlockID(0)
        problem_type = (DiscreteOrdinatesCurvilinearProblem if geometry == "sphere"
                        else DiscreteOrdinatesProblem)
        for method in ("classic_richardson", "petsc_richardson", "petsc_gmres"):
            problem = problem_type(
                mesh=grid,
                num_groups=1,
                groupsets=[
                    {
                        "groups_from_to": (0, 0),
                        "angular_quadrature": quad,
                        "inner_linear_method": method,
                        "l_abs_tol": 1.0e-14,
                        "l_max_its": 4,
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
                    print(f"CASE {geometry} {method} {op} rel_err {err:.6e}")

    if rank == 0:
        print(f"DSA_REFLECTING_RESIDUAL_MAX_REL_ERR {max_err:.6e}")
