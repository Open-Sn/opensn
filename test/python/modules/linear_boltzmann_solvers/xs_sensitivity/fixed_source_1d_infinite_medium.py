#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Semi-analytic 1D reflective infinite-medium sensitivity check."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.post import CrossSectionSensitivityPostprocessor
else:
    rank = 0
    barrier = MPIBarrier


def remove_rank_file(prefix):
    try:
        os.remove(prefix + str(rank) + ".h5")
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    sigma_t = 1.0
    scatter_ratio = 0.2
    sigma_s = sigma_t * scatter_ratio
    sigma_a = sigma_t - sigma_s
    q = 2.0
    length = 1.0
    num_cells = 24

    nodes = [i * length / num_cells for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, scatter_ratio)
    src = VolumetricSource(block_ids=[0], group_strength=[q])
    response_src = VolumetricSource(block_ids=[0], group_strength=[1.0])
    boundaries = [{"name": "zmin", "type": "reflecting"}, {"name": "zmax", "type": "reflecting"}]

    quad = GLProductQuadrature1DSlab(n_polar=32, scattering_order=0)
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-11,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=boundaries,
        options={
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )
    solver = SteadyStateSourceSolver(problem=phys)
    solver.Initialize()
    solver.Execute()

    fwd_phi_prefix = "xs_sens_inf_fwd_phi_"
    fwd_psi_prefix = "xs_sens_inf_fwd_psi_"
    phys.WriteFluxMoments(fwd_phi_prefix)
    phys.WriteAngularFluxes(fwd_psi_prefix)

    phys.SetAdjoint(True)
    phys.SetBoundaryOptions(boundary_conditions=boundaries)
    phys.SetVolumetricSources(volumetric_sources=[response_src])
    solver.Execute()
    adj_phi_prefix = "xs_sens_inf_adj_phi_"
    adj_psi_prefix = "xs_sens_inf_adj_psi_"
    phys.WriteFluxMoments(adj_phi_prefix)
    phys.WriteAngularFluxes(adj_psi_prefix)

    sigma_pp = CrossSectionSensitivityPostprocessor(
        problem=phys,
        sensitivity_type="sigma_t",
        forward_angular_fluxes=fwd_psi_prefix,
        adjoint_angular_fluxes=adj_psi_prefix,
    )
    sigma_pp.Execute()
    sigma_sens = sigma_pp.GetValue()[0][0]

    scatter_pp = CrossSectionSensitivityPostprocessor(
        problem=phys,
        sensitivity_type="scatter",
        moment=0,
        from_group=0,
        to_group=0,
        forward_flux_moments=fwd_phi_prefix,
        adjoint_flux_moments=adj_phi_prefix,
    )
    scatter_pp.Execute()
    scatter_sens = scatter_pp.GetValue()[0][0]

    analytic_magnitude = length * q / (sigma_a * sigma_a)
    sigma_err = abs(sigma_sens + analytic_magnitude)
    scatter_err = abs(scatter_sens - analytic_magnitude)

    if rank == 0:
        print(f"XS_SENS_INF_SIGMA_ERR={sigma_err:.12e}")
        print(f"XS_SENS_INF_SCATTER_ERR={scatter_err:.12e}")

    barrier()
    for prefix in [fwd_phi_prefix, fwd_psi_prefix, adj_phi_prefix, adj_psi_prefix]:
        remove_rank_file(prefix)
