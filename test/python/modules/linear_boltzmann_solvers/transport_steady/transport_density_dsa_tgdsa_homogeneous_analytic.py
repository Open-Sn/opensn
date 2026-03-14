#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2-group homogeneous analytic check with WGDSA and TGDSA enabled.

For uniform reflected media, group balance reduces to:
  A * phi = q,
  A[g,g'] = rho*sigma_t[g]*delta_{g,g'} - rho*S0[g'->g].
This test parses sigma_t and P0 transfer terms from XS, solves the 2x2 linear
system analytically, and compares to transport+DSA scalar fluxes.
The numeric group averages must match the 2x2 analytic solution.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def parse_xs_moment0(filename):
    sigma_t = {}
    s0 = {}
    in_sigma_t = False
    in_transfer = False

    with open(filename, "r", encoding="utf-8") as fin:
        for line in fin:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = raw.split()

            if parts[0] == "SIGMA_T_BEGIN":
                in_sigma_t = True
                continue
            if parts[0] == "SIGMA_T_END":
                in_sigma_t = False
                continue
            if parts[0] == "TRANSFER_MOMENTS_BEGIN":
                in_transfer = True
                continue
            if parts[0] == "TRANSFER_MOMENTS_END":
                in_transfer = False
                continue

            if in_sigma_t:
                g = int(parts[0])
                sigma_t[g] = float(parts[1])
            elif in_transfer and parts[0] == "M_GFROM_GTO_VAL":
                m = int(parts[1])
                if m != 0:
                    continue
                gfrom = int(parts[2])
                gto = int(parts[3])
                val = float(parts[4])
                s0[(gfrom, gto)] = val

    return [sigma_t[0], sigma_t[1]], s0


def solve_2x2(a00, a01, a10, a11, b0, b1):
    det = a00 * a11 - a01 * a10
    x0 = (b0 * a11 - a01 * b1) / det
    x1 = (a00 * b1 - b0 * a10) / det
    return x0, x1


def avg_phi(problem, group):
    ff = problem.GetScalarFluxFieldFunction()[group]
    lv = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("avg")
    ffi.SetLogicalVolume(lv)
    ffi.AddFieldFunction(ff)
    ffi.Initialize()
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/triangle_mesh2x2_cuts.obj",
        partitioner=KBAGraphPartitioner(nx=2, ny=2, nz=1, xcuts=[0.0], ycuts=[0.0]),
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)

    xs_file = "../transport_keigen/xs_water_g2.xs"
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_file)

    rho = 1.3
    q0, q1 = 1.0, 0.4

    pquad = GLCProductQuadrature2DXY(n_polar=6, n_azimuthal=24, scattering_order=1)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[{
            "groups_from_to": [0, 1],
            "angular_quadrature": pquad,
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-8,
            "l_max_its": 300,
            "apply_wgdsa": True,
            "apply_tgdsa": True,
            "wgdsa_l_abs_tol": 1.0e-8,
            "tgdsa_l_abs_tol": 1.0e-8,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[q0, q1])],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
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

    phi0_num = avg_phi(problem, 0)
    phi1_num = avg_phi(problem, 1)

    sigma_t, s0 = parse_xs_moment0(xs_file)

    # A * phi = q with A[g,g'] = rho*sigma_t[g]*delta - rho*S0[g'->g]
    a00 = rho * sigma_t[0] - rho * s0[(0, 0)]
    a01 = -rho * s0[(1, 0)]
    a10 = -rho * s0[(0, 1)]
    a11 = rho * sigma_t[1] - rho * s0[(1, 1)]

    phi0_ref, phi1_ref = solve_2x2(a00, a01, a10, a11, q0, q1)

    rel0 = abs(phi0_num - phi0_ref) / abs(phi0_ref)
    rel1 = abs(phi1_num - phi1_ref) / abs(phi1_ref)
    pass_flag = 1 if max(rel0, rel1) < 1.0e-2 else 0

    if rank == 0:
        print(f"DENSITY_DSA_REF_PHI0={phi0_ref:.8e}")
        print(f"DENSITY_DSA_REF_PHI1={phi1_ref:.8e}")
        print(f"DENSITY_DSA_NUM_PHI0={phi0_num:.8e}")
        print(f"DENSITY_DSA_NUM_PHI1={phi1_num:.8e}")
        print(f"DENSITY_DSA_REL0={rel0:.8e}")
        print(f"DENSITY_DSA_REL1={rel1:.8e}")
        print(f"DENSITY_DSA_ANALYTIC_PASS {pass_flag}")
