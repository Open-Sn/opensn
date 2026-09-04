#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Analytic 1D reflective infinite-medium k-eigenvalue and sensitivity check."""

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
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
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
    test_dir = os.path.dirname(__file__)
    sigma_t = 1.0
    sigma_s = 0.2
    # k = nu_sigma_f / (sigma_t - sigma_s) for a 1-group, homogeneous, infinite (fully reflected)
    # medium, so dk/d(nu_sigma_f) = 1 / (sigma_t - sigma_s) -- independent of the fissile XS
    # values (sigma_f=0.45, production=nu_sigma_f=0.9 in xs_sensitivity_fissile_1g.xs) themselves,
    # since k is linear in nu_sigma_f for a fixed medium.
    analytic_production_sens = 1.0 / (sigma_t - sigma_s)

    n_cells = 20
    nodes = [i / n_cells for i in range(n_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(os.path.join(test_dir, "xs_sensitivity_fissile_1g.xs"))

    quad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=0)
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
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )
    solver = PowerIterationKEigenSolver(problem=phys, max_iters=200, k_tol=1.0e-11)
    solver.Initialize()
    solver.Execute()
    k_eff = solver.GetEigenvalue()

    fwd_phi_prefix = "xs_sens_keigen_inf_fwd_phi_"
    adj_phi_prefix = "xs_sens_keigen_inf_adj_phi_"
    phys.WriteFluxMoments(fwd_phi_prefix)

    phys.SetAdjoint(True)
    adj_solver = PowerIterationKEigenSolver(problem=phys, max_iters=200, k_tol=1.0e-11)
    adj_solver.Initialize()
    adj_solver.Execute()
    adj_k_eff = adj_solver.GetEigenvalue()

    phys.WriteFluxMoments(adj_phi_prefix)

    production_pp = CrossSectionSensitivityPostprocessor(
        problem=phys,
        sensitivity_type="production",
        group=0,
        forward_flux_moments=fwd_phi_prefix,
        adjoint_flux_moments=adj_phi_prefix,
    )
    production_pp.Execute()
    production_pp.ApplyKEigenvalueScaling(k_eff)
    production_sens = production_pp.GetValue()[0][0]
    production_err = abs(production_sens - analytic_production_sens)

    if rank == 0:
        print(f"XS_SENS_KEIGEN_PRODUCTION_ERR={production_err:.12e}")

    barrier()
    for prefix in [fwd_phi_prefix, adj_phi_prefix]:
        remove_rank_file(prefix)
