#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Exact 3D check for all-group sigma_t angular-flux sensitivity."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.post import CrossSectionSensitivityPostprocessor
else:
    rank = 0
    barrier = MPIBarrier


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
    volume = 1.0

    meshgen = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(os.path.join(test_dir, "xs_sensitivity_2g.xs"))

    quad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=16, scattering_order=0)
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

    fwd_psi_prefix = "xs_sens_3d_fwd_psi_"
    fill_psi(phys, num_groups, [1.25, 2.5])
    phys.WriteAngularFluxes(fwd_psi_prefix)
    fill_psi(phys, num_groups, [3.0, 4.0])

    sigma_pp = CrossSectionSensitivityPostprocessor(
        problem=phys,
        sensitivity_type="sigma_t",
        forward_angular_fluxes=fwd_psi_prefix,
    )
    sigma_pp.Execute()
    sigma = sigma_pp.GetValue()[0]
    sigma_expected = [-volume * 1.25 * 3.0, -volume * 2.5 * 4.0]
    sigma_err = max(abs(a - b) for a, b in zip(sigma, sigma_expected))

    if rank == 0:
        print(f"XS_SENS_3D_SIGMA_ERR={sigma_err:.12e}")

    barrier()
    remove_rank_file(fwd_psi_prefix)
