#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Homogeneous reflected absorber/source scaling test.

For one-group pure absorption with reflecting boundaries:
  rho * sigma_t * phi = q  ->  phi = q / (rho * sigma_t).
Doubling the density halves scalar flux, so phi(rho=2)/phi(rho=1) = 0.5.
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


def max_phi(problem):
    ff = problem.GetScalarFluxFieldFunction()[0]
    lv = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("max")
    ffi.SetLogicalVolume(lv)
    ffi.AddFieldFunction(ff)
    ffi.Initialize()
    ffi.Execute()
    return ffi.GetValue()


def run_case(rho):
    nodes = [i / 20.0 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    sigma_t = 1.0
    sigma_s = 0.0
    xs.CreateSimpleOneGroup(sigma_t, sigma_s)

    src = VolumetricSource(block_ids=[0], group_strength=[1.0])
    quad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": quad,
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 200,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
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
    return max_phi(problem)


if __name__ == "__main__":
    phi1 = run_case(1.0)
    phi2 = run_case(2.0)

    ratio = phi2 / phi1 if phi1 != 0.0 else 0.0
    pass_flag = 1 if abs(ratio - 0.5) < 1.0e-4 else 0

    if rank == 0:
        print(f"DENSITY_REFLECT_RATIO={ratio:.8e}")
        print(f"DENSITY_REFLECT_PASS {pass_flag}")
