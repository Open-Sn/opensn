#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CSDA WGDSA is unaffected by an adjoint problem that shares the cross sections.

An adjoint problem puts the shared cross-section object in adjoint mode, which only changes the
getters that take no mode argument. The CSDA WGDSA operator must use the forward problem's mode,
so the forward solve is the same before and after an adjoint problem is built on the same object.

The problem is the 2D reflecting synthetic CSDA infinite medium of
transport_2d_cepxs_synthetic_csda_infinite_medium.py with WGDSA. With the cross sections' own
mode, the operator used the transposed in-groupset scattering after the adjoint problem was built,
and the solve took 24 instead of 17 sweeps. With no leakage the solution is spatially uniform, and
its group averages are the exact solution of the infinite-medium algebraic system of the discrete
CSDA model (compute_reference in transport_2d_cepxs_synthetic_csda_infinite_medium.py, solved here
in rational arithmetic). Checks:
  CSDA_WGDSA_SHARED_EXTRA_SWEEPS: sweeps after the adjoint problem minus sweeps before (0).
  CSDA_WGDSA_SHARED_MAX_PHI_REL_ERR: the largest relative error of the group averages of both
    solves against the exact solution (0 to the solver tolerance, classic Richardson with
    l_abs_tol 1e-10).
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


G = 3
# Exact infinite-medium group averages of the discrete CSDA model.
PHI_EXACT = [55 / 78, 3785 / 8073, 72011 / 1248624]


def group_averages(problem):
    averages = []
    for ff in problem.GetScalarFluxFieldFunction():
        ffi = FieldFunctionInterpolationVolume()
        ffi.SetOperationType("avg")
        ffi.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
        ffi.AddFieldFunction(ff)
        ffi.Execute()
        averages.append(ffi.GetValue())
    return averages


if __name__ == "__main__":
    xs = MultiGroupXS()
    xs_path = os.path.join(
        os.path.dirname(__file__), "../../../../assets/xs/cepxs_synthetic_csda_3g.bxslib"
    )
    xs.LoadFromCEPXS(
        xs_path,
        material_id=0,
        csda_format=True,
    )
    nodes = [i / 12 for i in range(13)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    quadrature = GLCProductQuadrature2DXY(n_polar=8, n_azimuthal=16, scattering_order=0)

    def forward_solve():
        problem = DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=G,
            groupsets=[
                {
                    "groups_from_to": (0, G - 1),
                    "angular_quadrature": quadrature,
                    "inner_linear_method": "classic_richardson",
                    "l_abs_tol": 1.0e-10,
                    "l_max_its": 2000,
                    "apply_wgdsa": True,
                    "wgdsa_l_abs_tol": 1.0e-12,
                }
            ],
            xs_map=[{"block_ids": [0], "xs": xs}],
            volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0, 0.0])],
            boundary_conditions=[
                {"name": name, "type": "reflecting"} for name in ("xmin", "xmax", "ymin", "ymax")
            ],
            options={"csda_enabled": True, "verbose_inner_iterations": False},
        )
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()
        return solver.GetNumSweeps(), group_averages(problem)

    sweeps_alone, phi_alone = forward_solve()
    DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=G,
        groupsets=[{"groups_from_to": (0, G - 1), "angular_quadrature": quadrature}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        options={"adjoint": True},
    )
    sweeps_shared, phi_shared = forward_solve()

    errors = [
        abs(phi[g] - PHI_EXACT[g]) / PHI_EXACT[g]
        for phi in (phi_alone, phi_shared)
        for g in range(G)
    ]
    # max() can skip a NaN, so any non-finite error is reported as infinity.
    max_rel_err = max(errors) if all(math.isfinite(e) for e in errors) else math.inf
    if rank == 0:
        print(f"CSDA_WGDSA_SHARED_EXTRA_SWEEPS={sweeps_shared - sweeps_alone}")
        print(f"CSDA_WGDSA_SHARED_MAX_PHI_REL_ERR={max_rel_err:.6e}")
