#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Disk-mesh density check with isotropic boundary.

A one-group pure absorber on a disk mesh with an isotropic boundary.
Attenuation should be stronger for larger rho: average scalar flux decreases
when density increases (phi(rho=2) < phi(rho=1)).
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
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


def run_case(rho):
    grid = FromFileMeshGenerator(filename="../../../../assets/mesh/disk.msh").Execute()
    grid.SetUniformBlockID(0)
    grid.SetUniformBoundaryID("outer")

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)

    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": pquad,
            "angle_aggregation_type": "single",
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 250,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[{"name": "outer", "type": "isotropic", "group_strength": [1.0]}],
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
    return avg_phi(problem)


if __name__ == "__main__":
    phi1 = run_case(1.0)
    phi2 = run_case(2.0)
    ratio = phi2 / phi1 if phi1 > 0.0 else 0.0
    pass_flag = 1 if (phi1 > 0.0 and phi2 > 0.0 and ratio < 1.0) else 0

    if rank == 0:
        print(f"DENSITY_DISK_PHI1={phi1:.8e}")
        print(f"DENSITY_DISK_PHI2={phi2:.8e}")
        print(f"DENSITY_DISK_RATIO={ratio:.8e}")
        print(f"DENSITY_DISK_PASS {pass_flag}")
