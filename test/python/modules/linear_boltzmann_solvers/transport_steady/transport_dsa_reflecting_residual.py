#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
WGDSA boundary residual on lagged reflecting boundaries.

An infinite medium is represented by an optically thin domain (0.1 mean free
paths across) with reflecting boundaries on every face. It holds a 1-group
material with scattering ratio c = 0.9999 and a uniform source. The total cross
section is 0.1 / width and the source is chosen so that the exact solution is
phi = Q / (sigma_t (1 - c)) = 1e4.

The option case_set selects the domains.
  "1d2d" (default): a 1D slab, a 2D square and a 1D sphere (reflecting outer
    surface), on orthogonal meshes 0.1 cm wide.
  "3d": a cube 0.1 cm wide on an orthogonal mesh, a 2 cm cube of unstructured
    tetrahedra (gmsh_all_tets.vtu), and a 2 cm cube of triangular prisms
    extruded from an unstructured triangle mesh (triangle_mesh2x2_cuts.obj).

The incoming angular flux on one boundary of each pair of opposing reflecting
boundaries (and on the reflecting surface of the sphere) is lagged by one
sweep. On such a boundary the error of a sweep has an incoming part, the change
of the lagged flux, which the WGDSA correction receives as an inward current.
Without it, the correction misses the particles that the lagged boundary has not
yet returned. In this thin domain that boundary term dominates.

Each case runs exactly max_its iterations (default 4), so the result measures
the acceleration, not the stopping test. With classic_richardson, the boundary
residual is the change of the lagged flux over a sweep. With petsc_richardson
and petsc_gmres, WGDSA is a preconditioner applied to a residual vector whose
lagged block is the boundary residual.

Measured errors after max_its iterations:
  1d2d, 2 ranks, 4 iterations: below 4e-6 with the residual. Without it the
    error is 23-59% with Richardson and 0.06-2.4% with GMRES.
  3d, 1 rank, 4 iterations: below 1.2e-5 with the residual. Without it the
    error is 65% with Richardson and 5e-4 to 6e-3 with GMRES.
  3d, 4 ranks, 8 iterations: below 6.2e-3 with the residual and 47% without
    it. On several ranks the sweep also lags angular fluxes on faces between
    partitions. WGDSA has no residual for those yet, and they limit the
    convergence of this case.

Check:
  DSA_REFLECTING_RESIDUAL_MAX_REL_ERR: max over the cases and inner solvers of
  |phi / 1e4 - 1| (min and max of phi over the domain).
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator, OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import (
        GLCProductQuadrature2DXY,
        GLCProductQuadrature3DXYZ,
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


def get_option(name, default):
    return globals()[name] if name in globals() else default


def orthogonal(dim, width, n, coord_sys="cartesian"):
    nodes = [i * width / n for i in range(n + 1)]
    return OrthogonalMeshGenerator(node_sets=[nodes] * dim, coord_sys=coord_sys).Execute()


def tetrahedra():
    return FromFileMeshGenerator(filename="../../../../assets/mesh/gmsh_all_tets.vtu").Execute()


def extruded_triangles():
    base = FromFileMeshGenerator(filename="../../../../assets/mesh/triangle_mesh2x2_cuts.obj")
    return ExtruderMeshGenerator(inputs=[base], layers=[{"z": 2.0, "n": 4}]).Execute()


if __name__ == "__main__":
    c, phi_exact = 0.9999, 1.0e4
    faces_1d = ("zmin", "zmax")
    faces_2d = ("xmin", "xmax", "ymin", "ymax")
    faces_3d = faces_2d + faces_1d

    def quad_1d():
        return GLProductQuadrature1DSlab(n_polar=8, scattering_order=0)

    def quad_2d():
        return GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8, scattering_order=0)

    def quad_3d():
        return GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

    # name, mesh, width of the domain (cm), quadrature, reflecting faces
    case_sets = {
        "1d2d": (
            ("slab", lambda: orthogonal(1, 0.1, 10), 0.1, quad_1d, faces_1d),
            ("xy", lambda: orthogonal(2, 0.1, 10), 0.1, quad_2d, faces_2d),
            ("sphere", lambda: orthogonal(1, 0.1, 10, "spherical"), 0.1,
             lambda: GLProductQuadrature1DSpherical(n_polar=8, scattering_order=0), ("zmax",)),
        ),
        "3d": (
            ("xyz", lambda: orthogonal(3, 0.1, 5), 0.1, quad_3d, faces_3d),
            ("tetrahedra", tetrahedra, 2.0, quad_3d, faces_3d),
            ("extruded_triangles", extruded_triangles, 2.0, quad_3d, faces_3d),
        ),
    }
    cases = case_sets[get_option("case_set", "1d2d")]
    max_its = int(get_option("max_its", 4))

    max_err = 0.0
    for geometry, make_grid, width, make_quad, names in cases:
        grid = make_grid()
        grid.SetUniformBlockID(0)
        if len(names) == 6:
            grid.SetOrthogonalBoundaries()
        problem_type = (DiscreteOrdinatesCurvilinearProblem if geometry == "sphere"
                        else DiscreteOrdinatesProblem)
        # The domain is 0.1 mean free paths wide, and phi = Q / (sigma_t (1 - c)) = phi_exact.
        sigma_t = 0.1 / width
        Q = phi_exact * sigma_t * (1.0 - c)
        xs = MultiGroupXS()
        xs.CreateSimpleOneGroup(sigma_t, c)
        quad = make_quad()
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
                        "l_max_its": max_its,
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
