#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SCDSA must not change the k-eigenvalue of a heterogeneous problem, for either diffusion
discretization.

A 2D two-group quarter core (fuel in [0, 20]^2 cm, water out to 30 cm, reflecting xmin/ymin,
vacuum xmax/ymax) is solved three times to tight tolerances: without acceleration, with SCDSA on
the PWLD (MIP) diffusion discretization, and with SCDSA on the continuous PWLC discretization.
Acceleration only changes how the power iteration converges, so all three eigenvalues must agree
to the iteration tolerance. The unaccelerated solve is the reference.

The PWLC variant restricts the discontinuous transport flux to a continuous field by nodal
averaging. Nodal averaging does not preserve the flux across the fuel/water interface, so the
diffusion solution must only correct the transport flux, never replace it.

Checks:
  SCDSA_PWLD_K_DIFF_PCM, SCDSA_PWLC_K_DIFF_PCM: (k_scdsa - k_reference) / k_reference * 1e5.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import (DiscreteOrdinatesProblem, PowerIterationKEigenSolver,
                                 SCDSAAcceleration)
    from pyopensn.logvol import RPPLogicalVolume


def solve(sdm):
    nodes = [30.0 * i / 24 for i in range(25)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=-1.0, xmax=20.0, ymin=-1.0, ymax=20.0, infz=True), 1, True)
    water = MultiGroupXS()
    water.LoadFromOpenSn("../../../../assets/xs/xs_water_g2.xs")
    fuel = MultiGroupXS()
    fuel.LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16,
                                                               scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 500,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": water}, {"block_ids": [1], "xs": fuel}],
        boundary_conditions=[{"name": "xmin", "type": "reflecting"},
                             {"name": "ymin", "type": "reflecting"}],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )
    kwargs = {}
    if sdm is not None:
        kwargs["acceleration"] = SCDSAAcceleration(problem=problem, sdm=sdm, l_abs_tol=1.0e-12,
                                                   pi_k_tol=1.0e-12)
    solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-11, max_iters=1000, **kwargs)
    solver.Initialize()
    solver.Execute()
    return solver.GetEigenvalue()


if __name__ == "__main__":
    k_reference = solve(None)
    k_pwld = solve("pwld")
    k_pwlc = solve("pwlc")
    if rank == 0:
        print(f"SCDSA_REFERENCE_K {k_reference:.10f}")
        print(f"SCDSA_PWLD_K_DIFF_PCM {(k_pwld - k_reference) / k_reference * 1.0e5:.6e}")
        print(f"SCDSA_PWLC_K_DIFF_PCM {(k_pwlc - k_reference) / k_reference * 1.0e5:.6e}")
