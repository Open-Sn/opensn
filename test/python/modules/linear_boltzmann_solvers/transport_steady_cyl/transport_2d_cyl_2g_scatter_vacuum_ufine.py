#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_2g_scatter_vacuum_ufine.py: 2-group scatter on gmsh fine
mesh with vacuum boundaries.
Expected: PASS=1 on comparison with orthogonal reference solution computed in
this input.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, OrthogonalMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume


if __name__ == "__main__":
    def run_problem(
        grid,
        xs,
        src,
        quad,
        angle_agg,
        inner_method,
        l_abs_tol,
        l_max_its,
        vol,
        boundary_conditions,
    ):
        problem = DiscreteOrdinatesCurvilinearProblem(
            mesh=grid,
            num_groups=2,
            groupsets=[
                {
                    "groups_from_to": (0, 1),
                    "angular_quadrature": quad,
                    "angle_aggregation_type": angle_agg,
                    "inner_linear_method": inner_method,
                    "l_abs_tol": l_abs_tol,
                    "l_max_its": l_max_its,
                }
            ],
            xs_map=[{"block_ids": [0], "xs": xs}],
            volumetric_sources=[src],
            boundary_conditions=boundary_conditions,
        )

        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()

        fflist = problem.GetScalarFluxFieldFunction(only_scalar_flux=False)

        ffi_g0 = FieldFunctionInterpolationVolume()
        ffi_g0.SetOperationType("max")
        ffi_g0.SetLogicalVolume(vol)
        ffi_g0.AddFieldFunction(fflist[0][0])
        ffi_g0.Initialize()
        ffi_g0.Execute()
        phi_g0 = ffi_g0.GetValue()

        ffi_g1 = FieldFunctionInterpolationVolume()
        ffi_g1.SetOperationType("max")
        ffi_g1.SetLogicalVolume(vol)
        ffi_g1.AddFieldFunction(fflist[1][0])
        ffi_g1.Initialize()
        ffi_g1.Execute()
        phi_g1 = ffi_g1.GetValue()

        return phi_g0, phi_g1

    level = "fine"
    mesh_file = "../../../../assets/mesh/rz_tri_fine.msh"
    ref_nr = 200
    ref_nz = 200
    tol = 5.0e-2

    meshgen = FromFileMeshGenerator(filename=mesh_file, coord_sys="cylindrical")
    grid = meshgen.Execute()
    vol = RPPLogicalVolume(xmin=0.0, xmax=1.0, ymin=0.0, ymax=2.0, infz=True)
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("transport_2d_cyl_2g_2.xs")
    src = VolumetricSource(block_ids=[0], group_strength=[1.5, 0.8])

    quad = GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=0)
    gmsh_bcs = [
        {"name": "rmin", "type": "reflecting"},
        {"name": "rmax", "type": "vacuum"},
        {"name": "zmin", "type": "vacuum"},
        {"name": "zmax", "type": "vacuum"},
    ]
    ref_bcs = [
        {"name": "rmin", "type": "reflecting"},
        {"name": "rmax", "type": "vacuum"},
        {"name": "zmin", "type": "vacuum"},
        {"name": "zmax", "type": "vacuum"},
    ]
    phi_g0, phi_g1 = run_problem(
        grid,
        xs,
        src,
        quad,
        angle_agg="single",
        inner_method="petsc_richardson",
        l_abs_tol=1.0e-12,
        l_max_its=5000,
        vol=vol,
        boundary_conditions=gmsh_bcs,
    )

    ref_nodes = [
        [i * (1.0 / ref_nr) for i in range(ref_nr + 1)],
        [i * (2.0 / ref_nz) for i in range(ref_nz + 1)],
    ]
    ref_meshgen = OrthogonalMeshGenerator(node_sets=ref_nodes, coord_sys="cylindrical")
    ref_grid = ref_meshgen.Execute()
    ref_grid.SetBlockIDFromLogicalVolume(vol, 0, True)
    phi_ref_g0, phi_ref_g1 = run_problem(
        ref_grid,
        xs,
        src,
        quad,
        angle_agg="single",
        inner_method="petsc_richardson",
        l_abs_tol=1.0e-12,
        l_max_its=200,
        vol=vol,
        boundary_conditions=ref_bcs,
    )

    ok_g0 = abs(phi_g0 - phi_ref_g0) / max(abs(phi_ref_g0), 1.0e-12) <= tol
    ok_g1 = abs(phi_g1 - phi_ref_g1) / max(abs(phi_ref_g1), 1.0e-12) <= tol
    ok = ok_g0 and ok_g1

    if rank == 0:
        print(f"G0_MAX {phi_g0:.12e}")
        print(f"G1_MAX {phi_g1:.12e}")
        print(f"G0_REF {phi_ref_g0:.12e}")
        print(f"G1_REF {phi_ref_g1:.12e}")
        print(f"TOL {tol:.6e}")
        print(f"PASS {int(ok)}")
