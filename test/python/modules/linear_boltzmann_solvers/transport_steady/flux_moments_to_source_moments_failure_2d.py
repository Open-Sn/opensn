#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A failed ReadFluxMomentsAndMakeSourceMoments leaves the flux vectors unchanged.

ReadFluxMomentsAndMakeSourceMoments(file_base, single_file) reads flux moments from file and makes
the external source moments from them. When making the source fails after the file has been read,
the call raises and must leave the problem's flux vectors (phi_new and phi_old) as they were.

An UncollidedProblem has no source function, so the source cannot be made for it: the call raises
after reading a compatible file. Two discrete-ordinates problems on the same mesh (one group, one
moment) with the uniform sources q = 1 and q = 3 write the flux files A and B. The UncollidedProblem
loads B into phi_old with ReadFluxMoments, so that phi_new (zero) and phi_old (B) differ from each
other and from A, and then reads A with ReadFluxMomentsAndMakeSourceMoments. The external source
moments of an UncollidedProblem are not read by anything, so they are not checked.

The error is a std::logic_error (a RuntimeError in Python) whose message names the problem and
MakeSourceMomentsFromPhi. The UncollidedProblem has its own name, so that the message is checked for
the problem's name and not for its type.

Checks:
  FLUX_MOMENTS_FAILED_CALL_FILE_DIFFERS: 1 if A differs from both phi_new and phi_old before the
    call, and phi_old is B (else the next checks are vacuous).
  FLUX_MOMENTS_FAILED_CALL_RAISED: 1 if the call raised a RuntimeError whose message names the
    problem and MakeSourceMomentsFromPhi.
  FLUX_MOMENTS_FAILED_CALL_PHI_NEW_UNCHANGED: 1 if phi_new is bit for bit the one before the call.
  FLUX_MOMENTS_FAILED_CALL_PHI_OLD_UNCHANGED: the same for phi_old.
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, UncollidedProblem
    from pyopensn.source import PointSource, VolumetricSource
    from pyopensn.xs import MultiGroupXS


def write_flux_file(grid, xs, q, file_base):
    """Solve a one-group problem with the uniform source q, write its flux moments, return phi."""
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=2, n_azimuthal=4, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[q])],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    problem.WriteFluxMoments(file_base)
    return np.array(problem.GetPhiNewLocal())


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    nodes = [0.25 * i for i in range(5)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.5)

    file_base_a = "flux_moments_to_source_moments_failure_2d_a_"
    file_base_b = "flux_moments_to_source_moments_failure_2d_b_"
    flux_a = write_flux_file(grid, xs, 1.0, file_base_a)
    flux_b = write_flux_file(grid, xs, 3.0, file_base_b)

    problem_name = "uncollided_point_source"
    uncollided = UncollidedProblem(
        name=problem_name,
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=[0.3, 0.4, 0.0], strength=[1.0])],
        near_source=[RPPLogicalVolume(xmin=-0.01, xmax=1.01, ymin=-0.01, ymax=1.01, infz=True)],
        scattering_order=0,
    )
    uncollided.ReadFluxMoments(file_base_b, False)

    phi_new = np.array(uncollided.GetPhiNewLocal())
    phi_old = np.array(uncollided.GetPhiOldLocal())
    file_differs = (
        flux_a.shape == phi_new.shape == phi_old.shape
        and phi_old.tobytes() == flux_b.tobytes()
        and bool(np.any(flux_a != phi_new))
        and bool(np.any(flux_a != phi_old))
    )

    raised = False
    try:
        uncollided.ReadFluxMomentsAndMakeSourceMoments(file_base_a, False)
    except Exception as error:
        message = str(error)
        raised = (
            type(error) is RuntimeError
            and problem_name in message
            and "MakeSourceMomentsFromPhi" in message
        )

    phi_new_unchanged = phi_new.tobytes() == np.array(uncollided.GetPhiNewLocal()).tobytes()
    phi_old_unchanged = phi_old.tobytes() == np.array(uncollided.GetPhiOldLocal()).tobytes()

    for file_base in (file_base_a, file_base_b):
        if os.path.exists(f"{file_base}{rank}.h5"):
            os.remove(f"{file_base}{rank}.h5")

    if rank == 0:
        print(f"FLUX_MOMENTS_FAILED_CALL_FILE_DIFFERS={int(file_differs)}")
        print(f"FLUX_MOMENTS_FAILED_CALL_RAISED={int(raised)}")
        print(f"FLUX_MOMENTS_FAILED_CALL_PHI_NEW_UNCHANGED={int(phi_new_unchanged)}")
        print(f"FLUX_MOMENTS_FAILED_CALL_PHI_OLD_UNCHANGED={int(phi_old_unchanged)}")
