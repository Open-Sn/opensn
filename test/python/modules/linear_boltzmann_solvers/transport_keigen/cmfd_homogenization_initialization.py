#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""CMFD fallback homogenization in voids and unequal-volume cells.

The Cartesian checks establish the documented bounded diffusion coefficient for
pure and near voids and the transport-rate result for a mixed void/material
coarse cell. The two-group check has identically zero transport flux in its second
group and fine-cell volumes in a 1:2 ratio; it therefore exercises volume fallback
and fails if fine cells are counted instead of integrated.

Each case runs one active power iteration with a negligible relaxation (1e-30),
so the final scalar flux is bitwise the transport iterate CMFD homogenized.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        CMFDAcceleration,
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
    )
    from pyopensn.xs import MultiGroupXS


def slab_groupset():
    return {
        "groups_from_to": (0, 0),
        "angular_quadrature": GLProductQuadrature1DSlab(n_polar=2, scattering_order=0),
        "inner_linear_method": "petsc_richardson",
        "l_max_its": 2,
        "l_abs_tol": 1.0e-12,
    }


def initialize_cartesian(materials, edges):
    grid = OrthogonalMeshGenerator(
        node_sets=[[float(i) for i in range(len(materials) + 1)]]
    ).Execute()
    grid.SetUniformBlockID(0)
    for block_id in range(1, len(materials)):
        grid.SetBlockIDFromLogicalVolume(
            RPPLogicalVolume(
                zmin=float(block_id), zmax=float(block_id + 1), infx=True, infy=True
            ),
            block_id,
            True,
        )
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[slab_groupset()],
        xs_map=[{"block_ids": [i], "xs": xs} for i, xs in enumerate(materials)],
        options={"use_precursors": False},
    )
    cmfd_options = {
        "problem": problem,
        "coarse_mesh": "identity",
        "relaxation": 1.0e-30,
    }
    if edges is not None:
        cmfd_options.update(
            coarse_mesh="local_aggregation",
            aggregation_size=len(materials),
            allow_mixed_material_aggregation=True,
        )
    cmfd = CMFDAcceleration(**cmfd_options)
    solver = PowerIterationKEigenSolver(
        problem=problem, acceleration=cmfd, max_iters=1, reset_phi0=True
    )
    solver.Initialize()
    solver.Execute()
    return problem, cmfd


if __name__ == "__main__":
    fissile = MultiGroupXS()
    fissile.LoadFromOpenSn("../../../../assets/xs/simple_fissile_1g.xs")
    void = MultiGroupXS()
    void.CreateSimpleOneGroup(sigma_t=0.0, c=0.0)
    near_void = MultiGroupXS()
    near_void.CreateSimpleOneGroup(sigma_t=1.0e-13, c=0.0)

    _, separated = initialize_cartesian([fissile, void, near_void], None)
    maximum_diffusion = 1.0 / (3.0e-12)
    void_error = max(
        abs(separated.GetHomogenizedDiffusionCoefficient(1)[0] - maximum_diffusion),
        abs(separated.GetHomogenizedDiffusionCoefficient(2)[0] - maximum_diffusion),
    ) / maximum_diffusion

    mixed_problem, mixed = initialize_cartesian([fissile, void], [0.0, 2.0])
    phi = mixed_problem.GetPhiNewLocal()
    material_flux = 0.5 * (phi[0] + phi[1])
    void_flux = 0.5 * (phi[2] + phi[3])
    expected_mixed_diffusion = (material_flux + void_flux) / (3.0 * material_flux)
    mixed_error = abs(
        mixed.GetHomogenizedDiffusionCoefficient(0)[0] - expected_mixed_diffusion
    )

    fallback_xs_file = "cmfd_homogenization_fallback_2g.xs"
    with open(fallback_xs_file, "w") as stream:
        stream.write(
            "NUM_GROUPS 2\nNUM_MOMENTS 1\n\n"
            "SIGMA_T_BEGIN\n0 1.0\n1 1.0\nSIGMA_T_END\n\n"
            "SIGMA_F_BEGIN\n0 0.5\n1 0.0\nSIGMA_F_END\n\n"
            "NU_PROMPT_BEGIN\n0 2.0\n1 0.0\nNU_PROMPT_END\n\n"
            "CHI_PROMPT_BEGIN\n0 1.0\n1 0.0\nCHI_PROMPT_END\n\n"
            "TRANSFER_MOMENTS_BEGIN\nTRANSFER_MOMENTS_END\n"
        )
    base = MultiGroupXS()
    base.LoadFromOpenSn(fallback_xs_file)
    dense = MultiGroupXS()
    dense.LoadFromOpenSn(fallback_xs_file)
    os.remove(fallback_xs_file)
    dense.Scale(2.0)
    fallback_grid = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0, 3.0]]).Execute()
    fallback_grid.SetUniformBlockID(0)
    fallback_grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(zmin=1.0, zmax=3.0, infx=True, infy=True), 1, True
    )
    fallback_problem = DiscreteOrdinatesProblem(
        mesh=fallback_grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLProductQuadrature1DSlab(
                    n_polar=2, scattering_order=0
                ),
                "inner_linear_method": "petsc_richardson",
                "l_max_its": 2,
                "l_abs_tol": 1.0e-12,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": base}, {"block_ids": [1], "xs": dense}],
        options={"use_precursors": False},
    )
    fallback_cmfd = CMFDAcceleration(
        problem=fallback_problem,
        coarse_mesh="local_aggregation",
        aggregation_size=2,
        allow_mixed_material_aggregation=True,
        relaxation=1.0e-30,
    )
    fallback_solver = PowerIterationKEigenSolver(
        problem=fallback_problem, acceleration=fallback_cmfd, max_iters=1
    )
    fallback_solver.Initialize()
    fallback_solver.Execute()

    expected_sigma_tr = (1.0 * 1.0 + 2.0 * 2.0) / 3.0
    expected_diffusion = 1.0 / (3.0 * expected_sigma_tr)
    fallback_errors = [
        abs(fallback_cmfd.GetHomogenizedRemovalXS(0)[1] - expected_sigma_tr),
        abs(fallback_cmfd.GetHomogenizedDiffusionCoefficient(0)[1] - expected_diffusion),
    ]
    fallback_error = max(fallback_errors)

    max_error = max(void_error, mixed_error, fallback_error)
    print(f"CMFD_HOMOGENIZATION_FALLBACK_MAX_ERROR {max_error:.16e}")
    if not math.isfinite(max_error) or max_error > 1.0e-12:
        raise RuntimeError(
            "CMFD initialization homogenization errors: "
            f"void={void_error}, mixed={mixed_error}, fallback={fallback_errors}"
        )
