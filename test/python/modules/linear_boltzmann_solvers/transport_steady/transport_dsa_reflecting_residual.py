#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
WGDSA residual of lagged angular fluxes: reflecting boundaries and broken sweep cycles.

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
  "cycles": meshes on which the sweep breaks cyclic dependencies whatever the
    partitioner library returns. A 0.1 cm square of 10 x 9 cells and a 0.1 cm
    cube of 5 x 5 x 5 cells are partitioned by cell number
    (LinearGraphPartitioner), which on 2 ranks cuts through a row of cells: the
    two partitions then depend on each other for half of the directions. A
    1 cm square of non-convex quadrilaterals (shestakov_quad_a0.05_n8.vtu) has
    cycles between the cells of one partition.
The option sphere_sweep_type selects "AAH" (default) or "CBC" for the sphere only.

The sweep lags an incoming angular flux by one sweep in two places: on one
boundary of each pair of opposing reflecting boundaries (and on the reflecting
surface of the sphere), and on the faces where it breaks a cyclic dependency,
between the cells of a partition or between partitions. There the error of a
sweep has an incoming part, the change of the lagged flux, which the WGDSA
correction receives as an inward current on the face of the cell that reads
the flux. Without it, the correction misses the particles that the lagged face
has not yet passed on. In this thin domain that term dominates.

Each case runs exactly max_its iterations (default 4), so the result measures
the acceleration, not the stopping test. With classic_richardson, the residual
is the change of the lagged flux over a sweep. With petsc_richardson and
petsc_gmres, WGDSA is a preconditioner applied to a residual vector whose
lagged block is that residual.

Measured errors after 4 iterations:
  1d2d, 2 ranks: below 4e-6 with the residuals. Without the residual of the
    reflecting boundaries the error is 23-59% with Richardson and 0.06-2.4%
    with GMRES. The 2D square has the same error whether the partitioner cuts
    it straight or not; without the residual of the broken cycles a cut with a
    step leaves 3.6%.
  3d, 1 rank: below 1.2e-5 with the residuals. Without the residual of the
    reflecting boundaries the error is 65% with Richardson and 5e-4 to 6e-3
    with GMRES.
  3d, 4 ranks: below 1.3e-5 with the residuals. Without the residual of the
    cycles broken between partitions the error is 1.5% to 3.5%.
  cycles, 2 ranks: below 1.2e-5 with the residuals. Without the residual of
    the broken cycles the error is 0.34% (square), 1.2% (cube) and 0.45%
    (non-convex quadrilaterals).

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
    from pyopensn.mesh import (
        ExtruderMeshGenerator,
        FromFileMeshGenerator,
        LinearGraphPartitioner,
        OrthogonalMeshGenerator,
    )
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


def cut_through_row(counts, width):
    """Orthogonal mesh partitioned by cell number, which cuts through a row of cells."""
    node_sets = [[i * width / n for i in range(n + 1)] for n in counts]
    return OrthogonalMeshGenerator(node_sets=node_sets,
                                   partitioner=LinearGraphPartitioner()).Execute()


def from_file(filename):
    return FromFileMeshGenerator(filename="../../../../assets/mesh/" + filename)


def tetrahedra():
    return from_file("gmsh_all_tets.vtu").Execute()


def extruded_triangles():
    base = from_file("triangle_mesh2x2_cuts.obj")
    return ExtruderMeshGenerator(inputs=[base], layers=[{"z": 2.0, "n": 4}]).Execute()


def nonconvex_quadrilaterals():
    grid = from_file("shestakov_quad_a0.05_n8.vtu").Execute()
    grid.SetOrthogonalBoundaries()
    return grid


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
        "cycles": (
            ("xy_cut_row", lambda: cut_through_row((10, 9), 0.1), 0.1, quad_2d, faces_2d),
            ("xyz_cut_row", lambda: cut_through_row((5, 5, 5), 0.1), 0.1, quad_3d, faces_3d),
            ("nonconvex_quadrilaterals", nonconvex_quadrilaterals, 1.0, quad_2d, faces_2d),
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
                sweep_type=(get_option("sphere_sweep_type", "AAH")
                            if geometry == "sphere" else "AAH"),
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
                **({"use_gpus": True} if get_option("use_gpus", False) else {}),
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
