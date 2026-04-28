#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Exact synthetic-buffer checks for 1D XS sensitivity kernels."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.post import CrossSectionSensitivityPostprocessor
else:
    rank = 0
    barrier = MPIBarrier


def fill_phi(problem, num_groups, num_moments, scale):
    phi = list(problem.GetPhiNewLocal())
    for i in range(len(phi)):
        g = i % num_groups
        m = (i // num_groups) % num_moments
        phi[i] = scale * (m + 1.0) * (g + 2.0)
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
    scattering_order = 1
    num_moments = 2
    length = 2.0

    meshgen = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0, length]])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(os.path.join(test_dir, "xs_sensitivity_2g.xs"))

    quad = GLProductQuadrature1DSlab(n_polar=8, scattering_order=scattering_order)
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
        xs_map=[{"block_ids": [0], "xs": xs}],
        options={"save_angular_flux": True},
    )
    solver = SteadyStateSourceSolver(problem=phys)
    solver.Initialize()

    fwd_phi_prefix = "xs_sens_1d_fwd_phi_"
    fwd_psi_prefix = "xs_sens_1d_fwd_psi_"
    fill_phi(phys, num_groups, num_moments, scale=1.0)
    fill_psi(phys, num_groups, [2.0, 3.0])
    phys.WriteFluxMoments(fwd_phi_prefix)
    phys.WriteAngularFluxes(fwd_psi_prefix)

    fill_phi(phys, num_groups, num_moments, scale=3.0)
    fill_psi(phys, num_groups, [5.0, 7.0])

    sigma_pp = CrossSectionSensitivityPostprocessor(
        problem=phys,
        sensitivity_type="sigma_t",
        forward_angular_fluxes=fwd_psi_prefix,
    )
    sigma_pp.Execute()
    sigma = sigma_pp.GetValue()[0]
    sigma_expected = [-length * 2.0 * 5.0, -length * 3.0 * 7.0]

    scatter_pp = CrossSectionSensitivityPostprocessor(
        problem=phys,
        sensitivity_type="scatter",
        from_group=0,
        to_group=1,
        forward_flux_moments=fwd_phi_prefix,
    )
    scatter_pp.Execute()
    scatter = scatter_pp.GetValue()[0]
    scatter_expected = [
        length * (3.0 * (0 + 1.0) * (1 + 2.0)) * ((0 + 1.0) * (0 + 2.0)),
        length * (3.0 * (1 + 1.0) * (1 + 2.0)) * ((1 + 1.0) * (0 + 2.0)),
    ]

    production_pp = CrossSectionSensitivityPostprocessor(
        problem=phys,
        sensitivity_type="production",
        from_group=0,
        to_group=1,
        forward_flux_moments=fwd_phi_prefix,
    )
    production_pp.Execute()
    production = production_pp.GetValue()[0][0]
    production_expected = scatter_expected[0]

    sigma_err = max(abs(a - b) for a, b in zip(sigma, sigma_expected))
    scatter_err = max(abs(a - b) for a, b in zip(scatter, scatter_expected))
    production_err = abs(production - production_expected)

    if rank == 0:
        print(f"XS_SENS_1D_SIGMA_ERR={sigma_err:.12e}")
        print(f"XS_SENS_1D_SCATTER_ERR={scatter_err:.12e}")
        print(f"XS_SENS_1D_PRODUCTION_ERR={production_err:.12e}")

    barrier()
    remove_rank_file(fwd_phi_prefix)
    remove_rank_file(fwd_psi_prefix)
