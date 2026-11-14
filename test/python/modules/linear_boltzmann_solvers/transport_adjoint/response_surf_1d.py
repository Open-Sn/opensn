#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Surface responses on a 1D slab, against volume duality (2 ranks).

1D quadratures pair opposite directions as (x, y, -z), which the adjoint reorientation and the
writer's adjoint boundary lookup use. The slab is 4 cm long with 40 cells, a source in z < 1, a
detector in z > 3, an interior plane at z = 2, a reflecting zmin, and P1 scattering. Two adjoint
problems are solved:
  A: driven by a volumetric source on the detector, with vacuum zmax.
  B: driven only by a direction- and group-dependent arbitrary incoming adjoint flux on zmax.
The reference is the volume duality <q, phi_adj> of the same adjoint data (EvaluateResponse).
  SURF_1D_A_PLANE_ERR, SURF_1D_B_PLANE_ERR: |R_plane / R_ref - 1|. The plane separates the source
    from the detector, so it carries the whole response (0 to the solver tolerance).
  SURF_1D_B_ZMAX_ERR: |R_zmax / R_ref - 1| for B, whose response enters through zmax.
  SURF_1D_A_ZMAX_FRACTION, SURF_1D_A_ZMIN_FRACTION: |R / R_ref| on the vacuum zmax and the
    reflecting zmin of A, which carry nothing (0).
A 1D orientation error would pair psi(mu) with psi_adj(-mu) and give errors of order 1.
"""


import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import AngularFluxFunction
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.response import ResponseEvaluator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier


if __name__ == "__main__":
    length, num_cells, plane = 4.0, 40, 2.0
    nodes = [i * length / num_cells for i in range(num_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(infx=True, infy=True, zmin=0.0, zmax=1.0), 1, True
    )
    detector = RPPLogicalVolume(infx=True, infy=True, zmin=3.0, zmax=length)
    grid.SetBlockIDFromLogicalVolume(detector, 2, True)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.6)
    xs_map = [{"block_ids": [0, 1, 2], "xs": xs}]
    quadrature = GLProductQuadrature1DSlab(n_polar=8, scattering_order=1)

    def make(adjoint, zmax):
        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[
                {
                    "groups_from_to": (0, 0),
                    "angular_quadrature": quadrature,
                    "inner_linear_method": "petsc_gmres",
                    "l_abs_tol": 1.0e-13,
                    "l_max_its": 1000,
                }
            ],
            xs_map=xs_map,
            boundary_conditions=[{"name": "zmin", "type": "reflecting"}, zmax],
            options={"save_angular_flux": True, "adjoint": adjoint},
        )

    def solve(problem):
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()

    vacuum = {"name": "zmax", "type": "vacuum"}
    interior = {"plane": ("z", plane)}
    boundaries = ["zmin", "zmax"]

    fwd = make(False, vacuum)
    fwd.SetVolumetricSources(
        volumetric_sources=[VolumetricSource(block_ids=[1], group_strength=[1.0])]
    )
    solve(fwd)
    fwd.WriteSurfaceAngularFluxes(
        "response_surf_1d_F_", boundary_surfaces=boundaries, interior_surfaces=interior
    )

    adj_a = make(True, vacuum)
    adj_a.SetVolumetricSources(
        volumetric_sources=[VolumetricSource(logical_volume=detector, group_strength=[0.4])]
    )
    solve(adj_a)
    adj_a.WriteFluxMoments("response_surf_1d_Aphi_")
    adj_a.WriteSurfaceAngularFluxes(
        "response_surf_1d_A_", boundary_surfaces=boundaries, interior_surfaces=interior
    )

    def incoming(group, direction):
        return 1.0 + 0.25 * direction + 0.3 * group

    adj_b = make(
        True, {"name": "zmax", "type": "arbitrary", "function": AngularFluxFunction(incoming)}
    )
    solve(adj_b)
    adj_b.WriteFluxMoments("response_surf_1d_Bphi_")
    adj_b.WriteSurfaceAngularFluxes(
        "response_surf_1d_B_", boundary_surfaces=boundaries, interior_surfaces=interior
    )

    surfaces = boundaries + ["plane"]

    def data(name, phi, surf):
        prefixes = {"surface_angular_fluxes": surf}
        if phi:
            prefixes["flux_moments"] = phi
        return {"name": name, "file_prefixes": prefixes, "surfaces": surfaces}

    ev = ResponseEvaluator(problem=fwd)
    ev.SetOptions(
        data=[
            data("F", None, "response_surf_1d_F_"),
            data("A", "response_surf_1d_Aphi_", "response_surf_1d_A_"),
            data("B", "response_surf_1d_Bphi_", "response_surf_1d_B_"),
        ],
        sources={"material": [{"block_id": 1, "strength": [1.0]}]},
    )
    qa = ev.EvaluateResponse("A")
    qb = ev.EvaluateResponse("B")
    ra = ev.EvaluateSurfaceResponse("F", "A", surfaces)
    rb = ev.EvaluateSurfaceResponse("F", "B", surfaces)

    if rank == 0:
        print(f"SURF_1D_A_PLANE_ERR={abs(ra['plane']['total'] / qa - 1.0):.6e}")
        print(f"SURF_1D_B_PLANE_ERR={abs(rb['plane']['total'] / qb - 1.0):.6e}")
        print(f"SURF_1D_B_ZMAX_ERR={abs(rb['zmax']['total'] / qb - 1.0):.6e}")
        print(f"SURF_1D_A_ZMAX_FRACTION={abs(ra['zmax']['total'] / qa):.6e}")
        print(f"SURF_1D_A_ZMIN_FRACTION={abs(ra['zmin']['total'] / qa):.6e}")

    barrier()
    if rank == 0:
        for name in os.listdir("."):
            if name.startswith("response_surf_1d_"):
                os.remove(name)
