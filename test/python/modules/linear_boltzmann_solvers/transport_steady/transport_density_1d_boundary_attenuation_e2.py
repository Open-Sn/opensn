#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Density attenuation check in a 1D pure absorber slab with isotropic inflow at zmin.

This test compares the simulated attenuation ratio at two depths to the semi-analytic
transport relation:
  phi(x) proportional to E2(rho * sigma_t * x),
  E2(tau) = integral_0^1 exp(-tau/mu) dmu.
Expected answer:
  phi(x2)/phi(x1) ~= E2(tau2)/E2(tau1), taui = rho * sigma_t * xi.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def avg_phi_over(problem, zmin, zmax):
    ff = problem.GetScalarFluxFieldFunction()[0]
    lv = RPPLogicalVolume(infx=True, infy=True, zmin=zmin, zmax=zmax)
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("avg")
    ffi.SetLogicalVolume(lv)
    ffi.AddFieldFunction(ff)
    ffi.Initialize()
    ffi.Execute()
    return ffi.GetValue()


def e2(tau):
    # Numerical integral for E2(tau) = integral_0^1 exp(-tau/mu) dmu
    n = 40000
    h = 1.0 / n
    s = 0.0
    for i in range(n):
        mu = (i + 0.5) * h
        s += math.exp(-tau / mu)
    return s * h


if __name__ == "__main__":
    nodes = [i / 160.0 for i in range(161)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    rho = 1.2
    sigma_t = 1.0

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, 0.0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": GLProductQuadrature1DSlab(n_polar=32, scattering_order=0),
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 250,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "zmin", "type": "isotropic", "group_strength": [1.0]},
            {"name": "zmax", "type": "vacuum"},
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

    x1 = 0.2
    x2 = 0.8
    dz = 0.03
    phi1 = avg_phi_over(problem, x1 - dz, x1 + dz)
    phi2 = avg_phi_over(problem, x2 - dz, x2 + dz)

    tau1 = rho * sigma_t * x1
    tau2 = rho * sigma_t * x2
    ratio_num = phi2 / phi1 if phi1 > 0.0 else 0.0
    ratio_ana = e2(tau2) / e2(tau1)
    rel = abs(ratio_num - ratio_ana) / ratio_ana if ratio_ana > 0.0 else 0.0

    pass_flag = 1 if rel < 5.0e-2 else 0

    if rank == 0:
        print(f"DENSITY_ATT_RATIO_NUM={ratio_num:.8e}")
        print(f"DENSITY_ATT_RATIO_ANA={ratio_ana:.8e}")
        print(f"DENSITY_ATT_REL={rel:.8e}")
        print(f"DENSITY_ATT_PASS {pass_flag}")
