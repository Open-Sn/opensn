#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Surface data hold the boundary values of the solve that produced the angular flux (2 ranks).

A one-group 2D problem is driven only by an isotropic incoming flux on xmin, with vacuum on the
other boundaries. An adjoint problem is driven by a volumetric source on the detector. By the
discrete surface identity, the `incoming` part of the xmin surface response is minus the forward
detector response QoI, and the plane x = 6, which separates xmin from the detector, carries QoI.

The forward surface data are written right after the solve, and again after SetBoundaryOptions
changes xmin to vacuum without a new solve. The file must pair the angular flux with the boundary
values of the solve, so both give the same responses. Taking xmin from the changed boundary gave an
incoming part of 0 instead of -QoI.
  SURF_BC_CHANGE_INCOMING_ERR: the largest |incoming(xmin) / (-QoI) - 1| of the two writes.
  SURF_BC_CHANGE_PLANE_ERR: the largest |R(plane) / QoI - 1| of the two writes.
All are 0 to the solver tolerance (1e-13).

With use_source_moments=True the sweep does not apply the incoming flux of non-reflecting
boundaries. A forward problem with the isotropic xmin source and source moments (the scattering
source of the first solve) therefore used no incoming flux, and its surface data must record none.
Recording the configured xmin flux gave an incoming part near -QoI.
  SURF_BC_CHANGE_MOMENTS_INCOMING: |incoming(xmin)| / QoI of that problem (exactly 0).

An angular flux that has not come from a completed solve has no boundary values to write, so the
write is rejected after ZeroPsi and after a mode change and back.
  SURF_BC_CHANGE_UNSOLVED_REJECTED: 1 when both writes are rejected and create no file on every
    rank.
With use_gpus=True the sweeps run on the GPU, which copies the angular flux back to the host before
the solver records the boundary values.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier

    def global_min(value):
        return MPI.COMM_WORLD.allreduce(value, op=MPI.MIN)

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.response import ResponseEvaluator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier

    def global_min(value):
        return MPIAllReduce(value, "min")

PREFIX = "response_surf_bc_change_"
USE_GPUS = globals().get("use_gpus", False)
SIGMA_DET = 0.3

if __name__ == "__main__":
    length, num_cells, plane = 10.0, 20, 6.0
    nodes = [i * length / num_cells for i in range(num_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    detector = RPPLogicalVolume(xmin=8.0, xmax=10.0, ymin=8.0, ymax=10.0, infz=True)
    grid.SetBlockIDFromLogicalVolume(detector, 2, True)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(0.3, 0.8)
    quadrature = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16, scattering_order=0)
    boundaries = ["xmin", "xmax", "ymin", "ymax"]
    interior = {"plane": ("x", plane)}

    def make_problem(adjoint, xmin, use_source_moments=False):
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
                    "gmres_restart_interval": 100,
                }
            ],
            xs_map=[{"block_ids": [0, 2], "xs": xs}],
            boundary_conditions=[xmin] + [{"name": b, "type": "vacuum"} for b in boundaries[1:]],
            options={
                "save_angular_flux": True,
                "adjoint": adjoint,
                "use_source_moments": use_source_moments,
            },
            **({"use_gpus": True} if USE_GPUS else {}),
        )

    def solve(problem):
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()

    def write(problem, name):
        problem.WriteSurfaceAngularFluxes(
            PREFIX + name, boundary_surfaces=boundaries, interior_surfaces=interior
        )

    def unsolved_rejected(problem, name):
        try:
            write(problem, name)
        except ValueError as error:
            created = os.path.exists(f"{PREFIX}{name}{rank}.h5")
            return "completed solve" in str(error) and not created
        return False

    forward = make_problem(False, {"name": "xmin", "type": "isotropic", "group_strength": [1.0]})
    solve(forward)
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("sum")
    ffi.SetLogicalVolume(detector)
    ffi.AddFieldFunction(forward.GetScalarFluxFieldFunction()[0])
    ffi.Execute()
    qoi = SIGMA_DET * ffi.GetValue()
    write(forward, "solved")
    forward.CreateAndWriteSourceMoments(PREFIX + "source_moments")
    forward.SetBoundaryOptions(boundary_conditions=[{"name": "xmin", "type": "vacuum"}])
    write(forward, "changed")

    forward.ZeroPsi()
    zeroed = unsolved_rejected(forward, "zeroed")
    solve(forward)
    forward.SetAdjoint(True)
    forward.SetAdjoint(False)
    switched = unsolved_rejected(forward, "switched")

    isotropic_xmin = {"name": "xmin", "type": "isotropic", "group_strength": [1.0]}
    moments = make_problem(False, isotropic_xmin, use_source_moments=True)
    moments.ReadSourceMoments(PREFIX + "source_moments", False)
    solve(moments)
    write(moments, "moments")

    adjoint = make_problem(True, {"name": "xmin", "type": "vacuum"})
    adjoint.SetVolumetricSources(
        volumetric_sources=[VolumetricSource(logical_volume=detector, group_strength=[SIGMA_DET])]
    )
    solve(adjoint)
    write(adjoint, "adjoint")

    surfaces = boundaries + ["plane"]
    evaluator = ResponseEvaluator(problem=adjoint)
    evaluator.SetOptions(
        data=[
            {
                "name": name,
                "file_prefixes": {"surface_angular_fluxes": PREFIX + name},
                "surfaces": surfaces,
            }
            for name in ("solved", "changed", "moments", "adjoint")
        ]
    )

    def finite_err(value):
        # max() can skip a NaN, so a non-finite error is reported as infinity.
        return value if math.isfinite(value) else math.inf

    incoming_err = 0.0
    plane_err = 0.0
    for name in ("solved", "changed"):
        responses = evaluator.EvaluateSurfaceResponse(name, "adjoint", surfaces)
        incoming = abs(-responses["xmin"]["incoming"] / qoi - 1.0)
        incoming_err = max(incoming_err, finite_err(incoming))
        plane_err = max(plane_err, finite_err(abs(responses["plane"]["total"] / qoi - 1.0)))

    responses = evaluator.EvaluateSurfaceResponse("moments", "adjoint", ["xmin"])
    moments_incoming = finite_err(abs(responses["xmin"]["incoming"] / qoi))

    # Each rank checks its own file, so every rank must reject both writes.
    unsolved_rejected_all = global_min(int(zeroed and switched))

    if rank == 0:
        print(f"SURF_BC_CHANGE_INCOMING_ERR={incoming_err:.6e}")
        print(f"SURF_BC_CHANGE_PLANE_ERR={plane_err:.6e}")
        print(f"SURF_BC_CHANGE_MOMENTS_INCOMING={moments_incoming:.6e}")
        print(f"SURF_BC_CHANGE_UNSOLVED_REJECTED={unsolved_rejected_all}")

    barrier()
    if rank == 0:
        for file_name in os.listdir("."):
            if file_name.startswith(PREFIX):
                os.remove(file_name)
