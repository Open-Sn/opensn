#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Leakage and flux trend check versus density for reflected and vacuum slabs.

This test verifies qualitative physics:
1. Reflecting boundaries: net leakage should be ~0.
2. Vacuum boundaries: increasing rho increases attenuation, reducing both leakage
   and interior scalar flux.
Expected answer:
  leak_reflect ~ 0 and leak_vac(rho=2) < leak_vac(rho=1), phi_vac(rho=2) < phi_vac(rho=1).
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


def run_case(rho, reflective):
    nodes = [i / 20.0 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)

    if reflective:
        bcs = [{"name": "zmin", "type": "reflecting"}, {"name": "zmax", "type": "reflecting"}]
    else:
        bcs = [{"name": "zmin", "type": "vacuum"}, {"name": "zmax", "type": "vacuum"}]

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 200,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=bcs,
        density={"default_density": rho},
        options={
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    leak = problem.ComputeLeakage(["zmin", "zmax"])
    total_leak = float(leak["zmin"][0] + leak["zmax"][0])
    mean_phi = sum(problem.GetPhiNewLocal()) / len(problem.GetPhiNewLocal())
    return total_leak, mean_phi


if __name__ == "__main__":
    leak_ref, phi_ref = run_case(1.0, True)
    leak_vac_1, phi_vac_1 = run_case(1.0, False)
    leak_vac_2, phi_vac_2 = run_case(2.0, False)

    leak_ratio = leak_vac_2 / leak_vac_1 if leak_vac_1 != 0.0 else 0.0
    phi_ratio = phi_vac_2 / phi_vac_1 if phi_vac_1 != 0.0 else 0.0

    pass_flag = 1 if (
        abs(leak_ref) < 1.0e-9
        and leak_vac_1 > 0.0
        and leak_vac_2 > 0.0
        and leak_vac_2 < leak_vac_1
        and phi_vac_2 < phi_vac_1
        and leak_ratio < 1.0
        and phi_ratio < 1.0
    ) else 0

    if rank == 0:
        print(f"DENSITY_LEAK_REFLECT={leak_ref:.8e}")
        print(f"DENSITY_LEAK_VAC_1={leak_vac_1:.8e}")
        print(f"DENSITY_LEAK_VAC_2={leak_vac_2:.8e}")
        print(f"DENSITY_LEAK_RATIO={leak_ratio:.8e}")
        print(f"DENSITY_LEAK_PHI_RATIO={phi_ratio:.8e}")
        print(f"DENSITY_LEAKAGE_PASS {pass_flag}")
