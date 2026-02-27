#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Synthetic 3-group CEPXS reader test with analytic infinite-medium solution.

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
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolation,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":
    # 1D slab mesh (aligned with z)
    num_cells = 50
    nodes = [i / num_cells for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Load synthetic CEPXS library
    xs = MultiGroupXS()
    xs.LoadFromCEPXS("cepxs_synthetic_3g.bxslib", material_id=0)

    # Volumetric sources (per group)
    q0 = 1.0
    q1 = 0.7
    q2 = 0.3
    mg_src = VolumetricSource(block_ids=[0], group_strength=[q0, q1, q2])

    # Angular quadrature
    pquad = GLProductQuadrature1DSlab(n_polar=32, scattering_order=3)

    # Physics
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=3,
        groupsets=[
            {
                "groups_from_to": (0, 2),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 200,
                "gmres_restart_interval": 50,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[mg_src],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={"field_function_prefix": "cepxs_synth3g"},
    )

    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Field functions (group-wise scalar flux)
    fflist = phys.GetScalarFluxFieldFunction(only_scalar_flux=False)

    # Average over entire domain
    lv = RPPLogicalVolume(xmin=-1.0, xmax=1.0, ymin=-1.0, ymax=1.0, zmin=-1.0, zmax=2.0)

    def avg_field(ff):
        ffi = FieldFunctionInterpolationVolume()
        ffi.SetOperationType("avg")
        ffi.SetLogicalVolume(lv)
        ffi.AddFieldFunction(ff)
        ffi.Initialize()
        ffi.Execute()
        return ffi.GetValue()

    phi0 = avg_field(fflist[0][0])
    phi1 = avg_field(fflist[1][0])
    phi2 = avg_field(fflist[2][0])

    # Analytic infinite-medium solution for the synthetic library
    # Sigma_t = [1.8, 1.4, 1.1]
    # Sigma_s (g_from -> g_to):
    #   0->0=0.6, 0->1=0.10, 0->2=0.04
    #   1->0=0.12, 1->1=0.5, 1->2=0.06
    #   2->0=0.05, 2->1=0.08, 2->2=0.4
    a11 = 1.8 - 0.6
    a12 = -0.12
    a13 = -0.05
    a21 = -0.10
    a22 = 1.4 - 0.5
    a23 = -0.08
    a31 = -0.04
    a32 = -0.06
    a33 = 1.1 - 0.4
    det = (
        a11 * (a22 * a33 - a23 * a32)
        - a12 * (a21 * a33 - a23 * a31)
        + a13 * (a21 * a32 - a22 * a31)
    )
    phi0_ref = (
        q0 * (a22 * a33 - a23 * a32)
        - a12 * (q1 * a33 - a23 * q2)
        + a13 * (q1 * a32 - a22 * q2)
    ) / det
    phi1_ref = (
        a11 * (q1 * a33 - a23 * q2)
        - q0 * (a21 * a33 - a23 * a31)
        + a13 * (a21 * q2 - q1 * a31)
    ) / det
    phi2_ref = (
        a11 * (a22 * q2 - q1 * a32)
        - a12 * (a21 * q2 - q1 * a31)
        + q0 * (a21 * a32 - a22 * a31)
    ) / det

    # Energy deposition (average over full domain)
    ff_edep = phys.CreateFieldFunction("energy_deposition", "energy_deposition")
    edep_avg = FieldFunctionInterpolationVolume()
    edep_avg.SetOperationType("avg")
    edep_avg.SetLogicalVolume(lv)
    edep_avg.AddFieldFunction(ff_edep)
    edep_avg.Initialize()
    edep_avg.Execute()
    edep = edep_avg.GetValue()

    # Reference energy deposition (sigma_edep dot phi_ref)
    edep_ref = 0.50 * phi0_ref + 0.30 * phi1_ref + 0.20 * phi2_ref

    if rank == 0:
        print(f"Phi0={phi0:.6e}")
        print(f"Phi1={phi1:.6e}")
        print(f"Phi2={phi2:.6e}")
        print(f"Edep={edep:.6e}")
        print(f"Phi0_ref={phi0_ref:.6e}")
        print(f"Phi1_ref={phi1_ref:.6e}")
        print(f"Phi2_ref={phi2_ref:.6e}")
        print(f"Edep_ref={edep_ref:.6e}")
