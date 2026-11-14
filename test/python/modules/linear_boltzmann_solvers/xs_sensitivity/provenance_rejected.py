#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The sensitivity postprocessor reads only forward files as forward and physical adjoint files as
adjoint.

The sensitivities pair forward and adjoint fluxes at the same direction. A file therefore has to
hold the kind of solution it is read as: a forward file has role FORWARD, and an adjoint file has
role ADJOINT and PHYSICAL orientation (written after the adjoint solver has reoriented the
solution). A file of another kind must raise an error instead of giving a wrong sensitivity.

The problem is the 1D reflecting infinite medium of fixed_source_1d_infinite_medium.py: one group,
sigma_t = 1, c = 0.2, source q = 2, length L = 1, and the response source 1. With
sigma_a = sigma_t (1 - c), the exact sensitivities are dR/dsigma_t = -L q / sigma_a^2 and
dR/dsigma_s = L q / sigma_a^2. The files read with the right roles must reproduce them.

Five kinds of file are written: forward phi and psi after the forward solve, adjoint phi and psi
after SetAdjoint(True) but before the adjoint solve (solver orientation), adjoint phi and psi after
the adjoint solve (physical orientation), and source moments (no solution role).

Checks:
  XS_SENS_PROVENANCE_SIGMA_ERR, XS_SENS_PROVENANCE_SCATTER_ERR: absolute errors of the sigma_t and
    scatter sensitivities from correct files.
  XS_SENS_PROVENANCE_SWAPPED_PSI_REJECTED: 1 if the forward and adjoint psi files read the wrong
    way round raise the role error.
  XS_SENS_PROVENANCE_SOLVER_PSI_REJECTED: 1 if a solver-oriented adjoint psi file raises the
    orientation error.
  XS_SENS_PROVENANCE_SWAPPED_PHI_REJECTED, XS_SENS_PROVENANCE_SOLVER_PHI_REJECTED: the same for
    flux moments.
  XS_SENS_PROVENANCE_SOURCE_MOMENTS_REJECTED: 1 if a source-moment file read as forward flux
    moments raises the role error.
  XS_SENS_PROVENANCE_CURRENT_STATE_REJECTED: 1 if leaving both prefixes empty, which would pair
    the current adjoint solution with itself, raises an error for the forward term.
"""

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
    barrier = MPIBarrier


def remove_rank_file(prefix):
    try:
        os.remove(prefix + str(rank) + ".h5")
    except FileNotFoundError:
        pass


def sensitivity(problem, sensitivity_type, forward, adjoint):
    """Return the sensitivity computed from the forward and adjoint file prefixes."""
    if sensitivity_type == "sigma_t":
        files = {"forward_angular_fluxes": forward, "adjoint_angular_fluxes": adjoint}
    else:
        files = {
            "moment": 0,
            "from_group": 0,
            "to_group": 0,
            "forward_flux_moments": forward,
            "adjoint_flux_moments": adjoint,
        }
    postprocessor = CrossSectionSensitivityPostprocessor(
        problem=problem, sensitivity_type=sensitivity_type, **files
    )
    postprocessor.Execute()
    return postprocessor.GetValue()[0][0]


def rejected(problem, sensitivity_type, forward, adjoint, expected, error_type=RuntimeError):
    """Return True if the sensitivity raises the error whose message contains expected."""
    try:
        sensitivity(problem, sensitivity_type, forward, adjoint)
    except error_type as error:
        return expected in str(error)
    except Exception:
        return False
    return False


if __name__ == "__main__":
    sigma_t = 1.0
    scatter_ratio = 0.2
    sigma_a = sigma_t * (1.0 - scatter_ratio)
    q = 2.0
    length = 1.0
    num_cells = 24

    nodes = [i * length / num_cells for i in range(num_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, scatter_ratio)
    src = VolumetricSource(block_ids=[0], group_strength=[q])
    response_src = VolumetricSource(block_ids=[0], group_strength=[1.0])
    boundaries = [{"name": "zmin", "type": "reflecting"}, {"name": "zmax", "type": "reflecting"}]

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=32, scattering_order=0),
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

    fwd_phi = "xs_sens_provenance_fwd_phi_"
    fwd_psi = "xs_sens_provenance_fwd_psi_"
    src_moments = "xs_sens_provenance_src_moments_"
    phys.WriteFluxMoments(fwd_phi)
    phys.WriteAngularFluxes(fwd_psi)
    phys.CreateAndWriteSourceMoments(src_moments)

    phys.SetAdjoint(True)
    phys.SetBoundaryOptions(boundary_conditions=boundaries)
    phys.SetVolumetricSources(volumetric_sources=[response_src])
    solver_phi = "xs_sens_provenance_solver_phi_"
    solver_psi = "xs_sens_provenance_solver_psi_"
    phys.WriteFluxMoments(solver_phi)
    phys.WriteAngularFluxes(solver_psi)

    solver.Execute()
    adj_phi = "xs_sens_provenance_adj_phi_"
    adj_psi = "xs_sens_provenance_adj_psi_"
    phys.WriteFluxMoments(adj_phi)
    phys.WriteAngularFluxes(adj_psi)

    analytic_magnitude = length * q / (sigma_a * sigma_a)
    sigma_err = abs(sensitivity(phys, "sigma_t", fwd_psi, adj_psi) + analytic_magnitude)
    scatter_err = abs(sensitivity(phys, "scatter", fwd_phi, adj_phi) - analytic_magnitude)

    role = "requested role"
    orientation = "requested orientation"
    swapped_psi = rejected(phys, "sigma_t", adj_psi, fwd_psi, role)
    solver_psi_rejected = rejected(phys, "sigma_t", fwd_psi, solver_psi, orientation)
    swapped_phi = rejected(phys, "scatter", adj_phi, fwd_phi, role)
    solver_phi_rejected = rejected(phys, "scatter", fwd_phi, solver_phi, orientation)
    source_moments = rejected(phys, "scatter", src_moments, adj_phi, role)
    current_state = rejected(phys, "scatter", "", "", "forward flux-moment file is required",
                             ValueError)

    if rank == 0:
        print(f"XS_SENS_PROVENANCE_SIGMA_ERR={sigma_err:.12e}")
        print(f"XS_SENS_PROVENANCE_SCATTER_ERR={scatter_err:.12e}")
        print(f"XS_SENS_PROVENANCE_SWAPPED_PSI_REJECTED={int(swapped_psi)}")
        print(f"XS_SENS_PROVENANCE_SOLVER_PSI_REJECTED={int(solver_psi_rejected)}")
        print(f"XS_SENS_PROVENANCE_SWAPPED_PHI_REJECTED={int(swapped_phi)}")
        print(f"XS_SENS_PROVENANCE_SOLVER_PHI_REJECTED={int(solver_phi_rejected)}")
        print(f"XS_SENS_PROVENANCE_SOURCE_MOMENTS_REJECTED={int(source_moments)}")
        print(f"XS_SENS_PROVENANCE_CURRENT_STATE_REJECTED={int(current_state)}")

    barrier()
    for prefix in [fwd_phi, fwd_psi, src_moments, solver_phi, solver_psi, adj_phi, adj_psi]:
        remove_rank_file(prefix)
