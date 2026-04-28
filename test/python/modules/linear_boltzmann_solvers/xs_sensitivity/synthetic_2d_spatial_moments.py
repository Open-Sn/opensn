#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Exact 2D checks for all-moment scattering and block-restricted sigma_t sensitivities."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.post import CrossSectionSensitivityPostprocessor
else:
    rank = 0
    barrier = MPIBarrier


def fill_phi(problem, num_groups, num_moments, scale):
    phi = list(problem.GetPhiNewLocal())
    for i in range(len(phi)):
        g = i % num_groups
        m = (i // num_groups) % num_moments
        phi[i] = scale * (m + 2.0) * (g + 1.0)
    problem.SetPhiNewLocal(phi)


def fill_psi(problem, num_groups, group_values):
    for psi_gs in problem.GetPsi():
        for i in range(len(psi_gs)):
            psi_gs[i] = group_values[i % num_groups]


def remove_rank_file(prefix):
    try:
        os.remove(prefix + str(rank) + ".h5")
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    test_dir = os.path.dirname(__file__)
    num_groups = 2
    scattering_order = 2
    num_moments = 6
    total_area = 2.0
    block_area = 1.0

    meshgen = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0, 2.0], [0.0, 1.0]])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    right_half = RPPLogicalVolume(xmin=1.0, xmax=2.0, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(right_half, 1, True)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(os.path.join(test_dir, "xs_sensitivity_2g.xs"))

    quad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16, scattering_order=scattering_order)
    moment_ells = [0, 1, 1, 2, 2, 2]
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": quad,
                "inner_linear_method": "petsc_gmres",
            }
        ],
        xs_map=[{"block_ids": [0, 1], "xs": xs}],
        options={"save_angular_flux": True},
    )
    solver = SteadyStateSourceSolver(problem=phys)
    solver.Initialize()

    fwd_phi_prefix = "xs_sens_2d_fwd_phi_"
    fwd_psi_prefix = "xs_sens_2d_fwd_psi_"
    fill_phi(phys, num_groups, num_moments, scale=1.0)
    fill_psi(phys, num_groups, [4.0, 6.0])
    phys.WriteFluxMoments(fwd_phi_prefix)
    phys.WriteAngularFluxes(fwd_psi_prefix)

    fill_phi(phys, num_groups, num_moments, scale=2.0)
    fill_psi(phys, num_groups, [0.5, 1.5])

    sigma_pp = CrossSectionSensitivityPostprocessor(
        problem=phys,
        sensitivity_type="sigma_t",
        group=1,
        block_ids=[1],
        forward_angular_fluxes=fwd_psi_prefix,
    )
    sigma_pp.Execute()
    sigma = sigma_pp.GetValue()[0][0]
    sigma_expected = -block_area * 6.0 * 1.5

    scatter_pp = CrossSectionSensitivityPostprocessor(
        problem=phys,
        sensitivity_type="scatter",
        from_group=0,
        to_group=1,
        forward_flux_moments=fwd_phi_prefix,
    )
    scatter_pp.Execute()
    scatter = scatter_pp.GetValue()[0]
    scatter_expected = []
    for ell in range(scattering_order + 1):
        total = 0.0
        for m, moment_ell in enumerate(moment_ells):
            if moment_ell == ell:
                fwd = (m + 2.0) * (0 + 1.0)
                adj = 2.0 * (m + 2.0) * (1 + 1.0)
                total += total_area * adj * fwd
        scatter_expected.append(total)

    sigma_err = abs(sigma - sigma_expected)
    scatter_err = max(abs(a - b) for a, b in zip(scatter, scatter_expected))

    if rank == 0:
        print(f"XS_SENS_2D_SIGMA_BLOCK_ERR={sigma_err:.12e}")
        print(f"XS_SENS_2D_SCATTER_MOMENTS_ERR={scatter_err:.12e}")

    barrier()
    remove_rank_file(fwd_phi_prefix)
    remove_rank_file(fwd_psi_prefix)
