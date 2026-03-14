#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2D unstructured reflected pure-absorber/source analytic-value check.

In a one-group uniformly reflected system:
  rho * sigma_t * phi = q.
The volume-average scalar flux phi = q / (rho * sigma_t).
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def avg_phi(problem, group=0):
    ff = problem.GetScalarFluxFieldFunction()[group]
    lv = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("avg")
    ffi.SetLogicalVolume(lv)
    ffi.AddFieldFunction(ff)
    ffi.Initialize()
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/triangle_mesh2x2_cuts.obj",
        partitioner=KBAGraphPartitioner(nx=2, ny=2, nz=1, xcuts=[0.0], ycuts=[0.0]),
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)

    sigma_t = 1.0
    c = 0.0
    rho = 1.7
    qval = 2.0

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, c)

    pquad = GLCProductQuadrature2DXY(n_polar=6, n_azimuthal=24, scattering_order=0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": pquad,
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 200,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[qval])],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
        ],
        density={"default_density": rho},
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    phi_avg = avg_phi(problem)
    phi_analytic = qval / (rho * sigma_t)
    rel = abs(phi_avg - phi_analytic) / phi_analytic
    pass_flag = 1 if rel < 5.0e-3 else 0

    if rank == 0:
        print(f"DENSITY_UNSTRUCT_ABS_PHI={phi_avg:.8e}")
        print(f"DENSITY_UNSTRUCT_ABS_ANALYTIC={phi_analytic:.8e}")
        print(f"DENSITY_UNSTRUCT_ABS_REL={rel:.8e}")
        print(f"DENSITY_UNSTRUCT_ABS_PASS {pass_flag}")
