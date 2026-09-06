#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Small 2D, 2-group fixed-source problem for visualizing scalar-flux and
reaction-rate-density field functions.

The domain is a 4 cm by 4 cm square split at x=2 cm:
  block 0: source/scattering material with g0 -> g1 downscatter
  block 1: pure absorber material

The PVTU output contains phi_g000_m00, phi_g001_m00, rr_total, rr_absorption,
rr_total_g0, rr_total_g1, rr_absorption_g0, and rr_absorption_g1.
"""

import math
import os
import sys
import glob

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.fieldfunc import FieldFunctionGridBased, FieldFunctionInterpolationVolume


def volume_sum(field_function, logical_volume):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("sum")
    ffi.SetLogicalVolume(logical_volume)
    ffi.AddFieldFunction(field_function)
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    length = 4.0
    n_cells = 40
    dx = length / n_cells
    nodes = [i * dx for i in range(n_cells + 1)]

    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()

    whole_domain = RPPLogicalVolume(xmin=0.0, xmax=length, ymin=0.0, ymax=length, infz=True)
    source_block = RPPLogicalVolume(xmin=0.0, xmax=2.0, ymin=0.0, ymax=length, infz=True)
    absorber_block = RPPLogicalVolume(xmin=2.0, xmax=length, ymin=0.0, ymax=length, infz=True)

    grid.SetBlockIDFromLogicalVolume(whole_domain, 0, True)
    grid.SetBlockIDFromLogicalVolume(absorber_block, 1, True)

    xs_source = MultiGroupXS()
    xs_source.LoadFromOpenSn("../../../../assets/xs/rr_2d_2g_source_mat.xs")

    xs_absorber = MultiGroupXS()
    xs_absorber.LoadFromOpenSn("../../../../assets/xs/rr_2d_2g_absorber_mat.xs")

    num_groups = 2
    source = VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0])
    quad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16, scattering_order=0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": quad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_source},
            {"block_ids": [1], "xs": xs_absorber},
        ],
        volumetric_sources=[source],
        boundary_conditions=[
            {"name": "xmin", "type": "vacuum"},
            {"name": "xmax", "type": "vacuum"},
            {"name": "ymin", "type": "vacuum"},
            {"name": "ymax", "type": "vacuum"},
        ],
        options={
            "verbose_inner_iterations": False,
        },
    )

    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    scalar_ffs = phys.GetScalarFluxFieldFunction(True)
    groupwise_rr_ffs = []

    rr_total_ff = phys.CreateFieldFunction("rr_total", "sigma_t")
    rr_abs_ff = phys.CreateFieldFunction("rr_absorption", "sigma_a")
    for g in range(num_groups):
        groupwise_rr_ffs.append(
            phys.CreateFieldFunction(f"rr_total_g{g}", "sigma_t", group=g)
        )
        groupwise_rr_ffs.append(
            phys.CreateFieldFunction(f"rr_absorption_g{g}", "sigma_a", group=g)
        )

    for rr_ff in [rr_total_ff, rr_abs_ff] + groupwise_rr_ffs:
        if rr_ff.CanUpdate():
            rr_ff.Update()

    output_base = "rr_2d_2g_fields"
    FieldFunctionGridBased.ExportMultipleToPVTU(
        [scalar_ffs[0], scalar_ffs[1], rr_total_ff, rr_abs_ff] + groupwise_rr_ffs,
        output_base,
    )

    sigma_t = {
        0: [1.0, 0.7],
        1: [1.5, 2.0],
    }
    sigma_a = {
        0: [0.3, 0.4],
        1: [1.5, 2.0],
    }

    manual_total = 0.0
    manual_abs = 0.0
    manual_group_total = [0.0 for _ in range(num_groups)]
    manual_group_abs = [0.0 for _ in range(num_groups)]
    for block_id, logical_volume in [(0, source_block), (1, absorber_block)]:
        for g in range(num_groups):
            phi_integral = volume_sum(scalar_ffs[g], logical_volume)
            total_rate = sigma_t[block_id][g] * phi_integral
            abs_rate = sigma_a[block_id][g] * phi_integral
            manual_total += total_rate
            manual_abs += abs_rate
            manual_group_total[g] += total_rate
            manual_group_abs[g] += abs_rate

    rr_total = volume_sum(rr_total_ff, whole_domain)
    rr_abs = volume_sum(rr_abs_ff, whole_domain)

    rel_err_total = abs((rr_total - manual_total) / manual_total)
    rel_err_abs = abs((rr_abs - manual_abs) / manual_abs)

    assert math.isclose(rr_total, manual_total, rel_tol=1.0e-10)
    assert math.isclose(rr_abs, manual_abs, rel_tol=1.0e-10)

    groupwise_rr_sums = {}
    for g in range(num_groups):
        rr_total_g = volume_sum(groupwise_rr_ffs[2 * g], whole_domain)
        rr_abs_g = volume_sum(groupwise_rr_ffs[2 * g + 1], whole_domain)
        rel_err_total_g = abs((rr_total_g - manual_group_total[g]) / manual_group_total[g])
        rel_err_abs_g = abs((rr_abs_g - manual_group_abs[g]) / manual_group_abs[g])
        assert math.isclose(rr_total_g, manual_group_total[g], rel_tol=1.0e-10)
        assert math.isclose(rr_abs_g, manual_group_abs[g], rel_tol=1.0e-10)
        groupwise_rr_sums[g] = {
            "total": rr_total_g,
            "abs": rr_abs_g,
            "rel_err_total": rel_err_total_g,
            "rel_err_abs": rel_err_abs_g,
        }

    if rank == 0:
        pvtu_path = output_base + ".pvtu"
        with open(pvtu_path, "r", encoding="utf-8") as f:
            pvtu_text = f.read()

        has_phi0 = 'Name="phi_g000_m00"' in pvtu_text
        has_phi1 = 'Name="phi_g001_m00"' in pvtu_text
        has_rr_total = 'Name="rr_total"' in pvtu_text
        has_rr_abs = 'Name="rr_absorption"' in pvtu_text
        has_groupwise_fields = True
        for g in range(num_groups):
            has_groupwise_fields = has_groupwise_fields and f'Name="rr_total_g{g}"' in pvtu_text
            has_groupwise_fields = (
                has_groupwise_fields and f'Name="rr_absorption_g{g}"' in pvtu_text
            )
        fields_exported = int(
            has_phi0 and has_phi1 and has_rr_total and has_rr_abs and has_groupwise_fields
        )

        print("Relative-error-total=", rel_err_total)
        print("Relative-error-abs=", rel_err_abs)
        for g in range(num_groups):
            print(f"Manual-total-rate-g{g}=", manual_group_total[g])
            print(f"Manual-abs-rate-g{g}=", manual_group_abs[g])
            print(f"Relative-error-total-g{g}=", groupwise_rr_sums[g]["rel_err_total"])
            print(f"Relative-error-abs-g{g}=", groupwise_rr_sums[g]["rel_err_abs"])
        print(f"RR2D2GFieldsExported={fields_exported}")

        for path in [pvtu_path] + glob.glob(output_base + "_*.vtu"):
            os.remove(path)
        print(f"Removed {output_base}.pvtu and {output_base}_*.vtu")
