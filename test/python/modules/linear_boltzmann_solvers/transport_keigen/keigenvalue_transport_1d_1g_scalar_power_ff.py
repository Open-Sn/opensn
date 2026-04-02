#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D 1G keigenvalue test that writes scalar flux, on-demand power, and a
generic XS-driven field function to the same PVTU/VTU output set.
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
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, NonLinearKEigenSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased

# Mesh
L = 100.0
n_cells = 50
dx = L / n_cells
nodes = [i * dx for i in range(n_cells + 1)]
meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
grid = meshgen.Execute()
grid.SetUniformBlockID(0)

# XS
xs_simple_fissile = MultiGroupXS()
xs_simple_fissile.LoadFromOpenSn("simple_fissile.xs")

# Problem
phys = DiscreteOrdinatesProblem(
    mesh=grid,
    num_groups=1,
    groupsets=[
        {
            "groups_from_to": (0, 0),
            "angular_quadrature": GLProductQuadrature1DSlab(n_polar=32, scattering_order=0),
            "inner_linear_method": "petsc_gmres",
            "l_max_its": 500,
            "l_abs_tol": 1.0e-8,
        }
    ],
    xs_map=[
        {
            "block_ids": [0],
            "xs": xs_simple_fissile,
        }
    ],
    options={
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": True,
        "power_default_kappa": 1.0,
    },
)

k_solver = NonLinearKEigenSolver(
    problem=phys,
    nl_max_its=5000,
    nl_abs_tol=1.0e-8,
)
k_solver.Initialize()
k_solver.Execute()

if not hasattr(phys, "GetScalarFluxFieldFunction"):
    raise RuntimeError("DiscreteOrdinatesProblem does not expose GetScalarFluxFieldFunction.")
if not hasattr(phys, "CreateFieldFunction"):
    raise RuntimeError("DiscreteOrdinatesProblem does not expose CreateFieldFunction.")

scalar_ff = phys.GetScalarFluxFieldFunction()[0]
power_ff = phys.CreateFieldFunction("power_generation", "power")
generic_ff = phys.CreateFieldFunction("fission_rate", "sigma_f", power_normalization_target=1.0)

out_base = "keigen_1d_scalar_power"
FieldFunctionGridBased.ExportMultipleToPVTU([scalar_ff, power_ff, generic_ff], out_base)

if rank == 0:
    pvtu_path = out_base + ".pvtu"
    with open(pvtu_path, "r", encoding="utf-8") as f:
        pvtu_text = f.read()
    n_arrays = pvtu_text.count("PDataArray")
    has_scalar = 'Name="phi_g000_m00"' in pvtu_text
    has_power = 'Name="power_generation"' in pvtu_text
    has_fission_rate = 'Name="fission_rate"' in pvtu_text
    same_file_ok = int(has_scalar and has_power and has_fission_rate)
    print(f"ScalarPowerAndXSSamePVTU={same_file_ok}")
    print(f"Wrote {pvtu_path} with {n_arrays} point-data arrays (expect >= 3).")
