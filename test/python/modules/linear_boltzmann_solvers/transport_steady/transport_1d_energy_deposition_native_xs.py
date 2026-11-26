#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""
Energy deposition read from an OpenSn cross-section file, in an infinite medium.

A 1D slab [0, 10] cm with 10 cells and reflecting faces holds a uniform volumetric source of
strength 1 in group 0. The two groups are a neutron group (sigma_t = 1, self-scatter 0.5) and a
photon group (sigma_t = 0.1, self-scatter 0.02) fed by the transfer 0 -> 1 of 0.3. The file gives
SIGMA_A (0.5, 0.08) and the ENERGY_DEPOSITION section (0.01, 0.002). With no leakage the flux is
uniform and exact:
  phi_0 = 1 / (1 - 0.5) = 2,
  phi_1 = 0.3 phi_0 / (0.1 - 0.02) = 7.5,
  energy deposition rate = 0.01 phi_0 + 0.002 phi_1 = 0.035 per cm^3,
so the integral of the "energy_deposition" field function over the slab is 0.35.

Checks: EDEP_INTEGRAL, EDEP_AVG, PHI_0_AVG and PHI_1_AVG, each against its exact value with a
tolerance of 1e-10 times that value.
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
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

XS_TEXT = """NUM_GROUPS 2
NUM_MOMENTS 1
SIGMA_T_BEGIN
0 1.0
1 0.1
SIGMA_T_END
SIGMA_A_BEGIN
0 0.5
1 0.08
SIGMA_A_END
TRANSFER_MOMENTS_BEGIN
M_GFROM_GTO_VAL 0 0 0 0.5
M_GFROM_GTO_VAL 0 0 1 0.3
M_GFROM_GTO_VAL 0 1 1 0.02
TRANSFER_MOMENTS_END
ENERGY_DEPOSITION_BEGIN
0 0.01
1 0.002
ENERGY_DEPOSITION_END
"""


def volume_value(ff, op_type, logical_volume):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType(op_type)
    ffi.SetLogicalVolume(logical_volume)
    ffi.AddFieldFunction(ff)
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    # Each rank reads its own copy of the file.
    xs_file = f"transport_1d_energy_deposition_native_xs_r{rank}.xs"
    with open(xs_file, "w") as f:
        f.write(XS_TEXT)
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_file)
    os.remove(xs_file)

    nodes = [float(i) for i in range(11)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 1000,
                "gmres_restart_interval": 30,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0])],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={"verbose_inner_iterations": False},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    vol = RPPLogicalVolume(infx=True, infy=True, infz=True)
    edep = problem.CreateFieldFunction("energy_deposition", "energy_deposition")
    edep_integral = volume_value(edep, "sum", vol)
    edep_avg = volume_value(edep, "avg", vol)
    phi = problem.GetScalarFluxFieldFunction(only_scalar_flux=True)
    phi_avg = [volume_value(phi[g], "avg", vol) for g in range(2)]

    if rank == 0:
        print(f"EDEP_INTEGRAL={edep_integral:.17e}")
        print(f"EDEP_AVG={edep_avg:.17e}")
        print(f"PHI_0_AVG={phi_avg[0]:.17e}")
        print(f"PHI_1_AVG={phi_avg[1]:.17e}")
