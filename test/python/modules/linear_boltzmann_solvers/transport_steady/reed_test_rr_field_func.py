#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
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
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume


def volume_sum(field_function, logical_volume):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("sum")
    ffi.SetLogicalVolume(logical_volume)
    ffi.AddFieldFunction(field_function)
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    widths = [2.0, 1.0, 2.0, 1.0, 2.0]
    nrefs = [200, 200, 200, 200, 200]
    nodes = [0.0]
    for width, nref in zip(widths, nrefs):
        dx = width / nref
        for _ in range(nref):
            nodes.append(nodes[-1] + dx)

    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()

    z_min = 0.0
    for mat_id, width in enumerate(widths):
        z_max = z_min + width
        lv = RPPLogicalVolume(infx=True, infy=True, zmin=z_min, zmax=z_max)
        grid.SetBlockIDFromLogicalVolume(lv, mat_id, True)
        z_min = z_max

    sigt = [50.0, 5.0, 0.0, 1.0, 1.0]
    c = [0.0, 0.0, 0.0, 0.9, 0.9]
    xs_map = []
    for mat_id in range(len(widths)):
        xs = MultiGroupXS()
        xs.CreateSimpleOneGroup(sigt[mat_id], c[mat_id])
        xs_map.append({"block_ids": [mat_id], "xs": xs})

    src0 = VolumetricSource(block_ids=[0], group_strength=[50.0])
    src1 = VolumetricSource(block_ids=[3], group_strength=[1.0])
    gl_quad = GLProductQuadrature1DSlab(n_polar=128, scattering_order=0)

    num_groups = 1
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": gl_quad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-9,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=xs_map,
        volumetric_sources=[src0, src1],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
    )

    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    logvol_whole_domain = RPPLogicalVolume(infx=True, infy=True, zmin=0.0, zmax=sum(widths))

    total_rr_ff = phys.CreateFieldFunction("rr_total", "sigma_t")
    abs_rr_ff = phys.CreateFieldFunction("rr_abs", "sigma_a")
    block3_total_rr_ff = phys.CreateFieldFunction("rr_total_block3", "sigma_t", block_ids=[3])
    phi_ff = phys.GetScalarFluxFieldFunction(True)[0]

    rr_total = volume_sum(total_rr_ff, logvol_whole_domain)
    rr_abs = volume_sum(abs_rr_ff, logvol_whole_domain)
    rr_total_block3 = volume_sum(block3_total_rr_ff, logvol_whole_domain)

    manual_total = 0.0
    manual_abs = 0.0
    manual_total_block3 = 0.0
    for mat_id, width in enumerate(widths):
        zmin = float(sum(widths[:mat_id]))
        zmax = zmin + width
        lv = RPPLogicalVolume(infx=True, infy=True, zmin=zmin, zmax=zmax)
        phi_integral = volume_sum(phi_ff, lv)

        manual_total += sigt[mat_id] * phi_integral
        manual_abs += sigt[mat_id] * (1.0 - c[mat_id]) * phi_integral
        if mat_id == 3:
            manual_total_block3 += sigt[mat_id] * phi_integral

    rel_err_total = abs((rr_total - manual_total) / manual_total)
    rel_err_abs = abs((rr_abs - manual_abs) / manual_abs)
    rel_err_total_block3 = abs((rr_total_block3 - manual_total_block3) / manual_total_block3)

    assert math.isclose(rr_total, manual_total, rel_tol=1.0e-10)
    assert math.isclose(rr_abs, manual_abs, rel_tol=1.0e-10)
    assert math.isclose(rr_total_block3, manual_total_block3, rel_tol=1.0e-10)
    assert total_rr_ff.CanUpdate()
    total_rr_ff.Update()

    if rank == 0:
        print("Relative-error-total=", rel_err_total)
        print("Relative-error-abs=", rel_err_abs)
        print("Relative-error-total-block3=", rel_err_total_block3)
