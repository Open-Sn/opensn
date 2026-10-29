#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D slab k-eigenvalue solve using OpenMC-generated delayed-neutron MGXS data.

slab_delayed_mgxs.h5 (two groups, fuel with six delayed groups, water) was generated for this
geometry: fuel on [-2, 2], water on [2, 10], vacuum at both ends. In OpenMC the transverse faces
are reflective, so the problem is effectively 1D.

References for k:
- OpenMC multi-group mode on the same file: 0.35908 +/- 0.00002 (8.4e8 active histories; the
  collision, track-length and absorption estimators agree).
- OpenSn converged in space and angle (3840 cells, S512): 0.3586326.
- This test (120 cells, S8): 0.3581939. S8 is 4.4e-4 below the converged value in this leaky
  slab (leakage fraction 0.66), so the gold is a regression baseline.

OpenSn and OpenMC differ by 4.4e-4 because they read the file's absorption differently. OpenSn
recomputes absorption as total minus the nu-scatter row sum. OpenMC uses the file's absorption
dataset, which differs from that by statistical noise in the generated library (fast water
0.00040 against 0.00052). Scattering built from the file's absorption the way OpenMC does gives
an OpenSn k of 0.3590775, within 3e-6 of OpenMC. The multiplicity matrix changes k by about 1e-6,
since scatter_matrix is already nu-scatter.

The script prints:
- imported delayed-data summaries for each material
- OpenSn k_eff
- a global balance table
- per-group scalar-flux maxima in the fuel and water regions
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver


def xs_summary(name, xs):
    nu_sigma_f = list(xs.nu_sigma_f) if xs.is_fissionable else []
    nu_prompt_sigma_f = list(xs.nu_prompt_sigma_f) if xs.is_fissionable else []
    nu_delayed_sigma_f = list(xs.nu_delayed_sigma_f) if xs.is_fissionable else []
    return {
        "dataset": name,
        "is_fissionable": xs.is_fissionable,
        "num_groups": xs.num_groups,
        "num_precursors": xs.num_precursors,
        "sum_nu_sigma_f": sum(nu_sigma_f),
        "sum_nu_prompt_sigma_f": sum(nu_prompt_sigma_f),
        "sum_nu_delayed_sigma_f": sum(nu_delayed_sigma_f),
    }


def region_flux_max(field_function, logical_volume):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("max")
    ffi.SetLogicalVolume(logical_volume)
    ffi.AddFieldFunction(field_function)
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    xs_fuel = MultiGroupXS()
    xs_fuel.LoadFromOpenMC("slab_delayed_mgxs.h5", "fuel", 294.0)

    xs_water = MultiGroupXS()
    xs_water.LoadFromOpenMC("slab_delayed_mgxs.h5", "water", 294.0)

    if xs_fuel.num_groups != xs_water.num_groups:
        raise RuntimeError("Fuel and water MGXS datasets must have the same number of groups.")

    slab_min = -2.0
    fuel_max = 2.0
    slab_max = 10.0
    n_cells = 120
    dz = (slab_max - slab_min) / n_cells
    nodes = [slab_min + i * dz for i in range(n_cells + 1)]

    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(1)

    fuel_region = RPPLogicalVolume(infx=True, infy=True, zmin=slab_min, zmax=fuel_max)
    water_region = RPPLogicalVolume(infx=True, infy=True, zmin=fuel_max, zmax=slab_max)
    grid.SetBlockIDFromLogicalVolume(fuel_region, 0, True)

    num_groups = xs_fuel.num_groups
    pquad = GLProductQuadrature1DSlab(
        n_polar=8,
        scattering_order=max(xs_fuel.scattering_order, xs_water.scattering_order),
    )

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 500,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_fuel},
            {"block_ids": [1], "xs": xs_water},
        ],
        boundary_conditions=[
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={
            "use_precursors": True,
            "verbose_inner_iterations": True,
            "verbose_outer_iterations": True,
        },
    )

    solver = PowerIterationKEigenSolver(
        problem=problem,
        max_iters=200,
        k_tol=1.0e-8,
    )
    solver.Initialize()
    solver.Execute()

    balance = solver.ComputeBalanceTable()
    scalar_flux = problem.GetScalarFluxFieldFunction()

    if rank == 0:
        print("FUEL_XS_SUMMARY", xs_summary("fuel", xs_fuel))
        print("WATER_XS_SUMMARY", xs_summary("water", xs_water))
        print(f"OPENSN_KEFF {solver.GetEigenvalue():.8f}")
        print(f"OPENSN_BALANCE {dict(balance)}")

    for g in range(num_groups):
        fuel_phi_max = region_flux_max(scalar_flux[g], fuel_region)
        water_phi_max = region_flux_max(scalar_flux[g], water_region)
        if rank == 0:
            print(f"GROUP_{g}_FUEL_PHI_MAX {fuel_phi_max:.10e}")
            print(f"GROUP_{g}_WATER_PHI_MAX {water_phi_max:.10e}")
