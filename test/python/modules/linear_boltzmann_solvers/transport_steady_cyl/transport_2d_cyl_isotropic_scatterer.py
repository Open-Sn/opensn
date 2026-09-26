#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_isotropic_scatterer.py: heterogeneous RZ pure scatterer enclosed by isotropic
boundaries.

An inner region with sigma_t = 5 per cm (10 mean free paths across) sits in an outer region with
sigma_t = 1 per cm. Both are pure scatterers (sigma_s0 = sigma_t) with a P1 scattering moment of
0.3 sigma_t. The axis at r = 0 is reflecting and rmax, zmin, and zmax carry isotropic incoming
angular flux X.

psi = X in every direction is the exact solution: the flux is isotropic, so the P1 moment and
the angular-redistribution term vanish, and sigma_t X equals the scattering source sigma_s0 phi
when phi = X. OpenSn normalizes quadrature weights to sum to one, so phi = X everywhere. The
result depends on the boundary input, the scattering source, and the curvilinear redistribution
term using the same normalization; with the former 2 pi weight sum, phi would be 2 pi X.
Expected: PHI_MIN = PHI_MAX = 1.75.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
else:
    barrier = MPIBarrier


def write_pure_scatterer_xs(path, sigma_t):
    with open(path, "w") as f:
        f.write(f"NUM_GROUPS 1\nNUM_MOMENTS 2\n\nSIGMA_T_BEGIN\n0 {sigma_t!r}\nSIGMA_T_END\n\n")
        f.write(
            f"TRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 {sigma_t!r}\n"
            f"M_GFROM_GTO_VAL 1 0 0 {0.3 * sigma_t!r}\nTRANSFER_MOMENTS_END\n"
        )


if __name__ == "__main__":
    boundary_flux = 1.75
    nodes = [i * 0.25 for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes], coord_sys="cylindrical").Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.0, xmax=1.0, ymin=0.5, ymax=1.5, infz=True), 1, True
    )

    xs_files = ["transport_2d_cyl_isotropic_scatterer_outer.xs",
                "transport_2d_cyl_isotropic_scatterer_inner.xs"]
    if rank == 0:
        write_pure_scatterer_xs(xs_files[0], 1.0)
        write_pure_scatterer_xs(xs_files[1], 5.0)
    barrier()
    xs = []
    for path in xs_files:
        x = MultiGroupXS()
        x.LoadFromOpenSn(path)
        xs.append(x)
    barrier()
    if rank == 0:
        for path in xs_files:
            os.remove(path)

    quad = GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=1)
    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
                "gmres_restart_interval": 100,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs[0]}, {"block_ids": [1], "xs": xs[1]}],
        boundary_conditions=[{"name": "rmin", "type": "reflecting"}]
        + [
            {"name": name, "type": "isotropic", "group_strength": [boundary_flux]}
            for name in ("rmax", "zmin", "zmax")
        ],
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    everywhere = RPPLogicalVolume(infx=True, infy=True, infz=True)
    values = {}
    for op in ("min", "max"):
        ffi = FieldFunctionInterpolationVolume()
        ffi.SetOperationType(op)
        ffi.SetLogicalVolume(everywhere)
        ffi.AddFieldFunction(problem.GetScalarFluxFieldFunction()[0])
        ffi.Execute()
        values[op] = ffi.GetValue()

    if rank == 0:
        print(f"PHI_MIN={values['min']:.12e}")
        print(f"PHI_MAX={values['max']:.12e}")
