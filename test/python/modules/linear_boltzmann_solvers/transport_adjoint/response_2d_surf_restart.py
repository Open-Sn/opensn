#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D surface response from adjoint surface data written after a restart.

The model is that of response_2d_surf.py: the interior plane x = 6 separates the source from the
detector and all exterior boundaries are vacuum, so the surface response across the plane,

    QoI = int_S int_4pi (Omega . n) psi(r, Omega) psi^dagger(r, Omega) dOmega dA,

equals the volumetric forward QoI exactly for the upwind DFEM discretization.

An adjoint solve leaves the adjoint angular flux in the physical orientation, while restart data
are loaded in the internal sweep orientation, which stores psi^dagger(-Omega) under Omega. Surface
data written in that state would pair each forward direction with the opposite adjoint direction.
This test checks that

* the surface response from an adjoint solve that writes restart data equals the forward QoI;
* writing adjoint surface data right after a restart read is rejected without creating a file;
* executing the adjoint solver after the restart read, then writing, reproduces the forward QoI.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import KBAGraphPartitioner, OrthogonalMeshGenerator
    from pyopensn.response import ResponseEvaluator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier


if __name__ == "__main__":
    num_procs = 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Create mesh and assign material IDs.
    num_cells = 20
    length = 10.0
    cell_size = length / num_cells
    nodes = [i * cell_size for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes],
        partitioner=KBAGraphPartitioner(nx=1, ny=2, nz=1, ycuts=[length / 2.0]),
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)

    source_volume = RPPLogicalVolume(
        xmin=nodes[0], xmax=nodes[2], ymin=nodes[0], ymax=nodes[2], infz=True
    )
    grid.SetBlockIDFromLogicalVolume(source_volume, 1, True)

    detector_volume = RPPLogicalVolume(
        xmin=nodes[-2], xmax=nodes[-1], ymin=nodes[-2], ymax=nodes[-1], infz=True
    )
    grid.SetBlockIDFromLogicalVolume(detector_volume, 2, True)

    # Define cross sections.
    xs_background = MultiGroupXS()
    xs_background.CreateSimpleOneGroup(0.1, 0.9)

    xs_source = MultiGroupXS()
    xs_source.CreateSimpleOneGroup(2.0, 1.0)

    xs_detector = MultiGroupXS()
    xs_detector.CreateSimpleOneGroup(0.8, 0.0)

    xs_map = [
        {"block_ids": [0], "xs": xs_background},
        {"block_ids": [1], "xs": xs_source},
        {"block_ids": [2], "xs": xs_detector},
    ]
    source_block_id = 1
    detector_sigma_a = 0.8
    source_strength = 1.0
    source_area = (nodes[2] - nodes[0]) ** 2

    quadrature = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=32, scattering_order=0)

    def make_problem(options):
        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[
                {
                    "groups_from_to": (0, 0),
                    "angular_quadrature": quadrature,
                    "inner_linear_method": "petsc_gmres",
                    "l_abs_tol": 1.0e-12,
                    "l_max_its": 500,
                    "gmres_restart_interval": 100,
                },
            ],
            xs_map=xs_map,
            boundary_conditions=[
                {"name": name, "type": "vacuum"} for name in ("xmin", "xmax", "ymin", "ymax")
            ],
            options={"save_angular_flux": True, **options},
        )

    def make_adjoint_problem(options):
        problem = make_problem({"adjoint": True, **options})
        problem.SetVolumetricSources(
            volumetric_sources=[
                VolumetricSource(logical_volume=detector_volume, group_strength=[detector_sigma_a])
            ]
        )
        return problem

    restart_stem = "response_2d_surf_restart_"
    forward_surface_prefix = "RestartFwdSurf_p"
    adjoint_surface_prefix = "RestartAdjSurf_p"
    rejected_surface_prefix = "RestartRejectedSurf_p"
    restarted_surface_prefix = "RestartedAdjSurf_p"
    interior_surfaces = {"inter_x": ("x", 6.0)}
    # The interior surface name reads both stored sides, inter_x_u and inter_x_d.
    interior_tags = ["inter_x"]

    # Forward solve, volumetric QoI, and forward surface data.
    forward_problem = make_problem({})
    forward_problem.SetVolumetricSources(
        volumetric_sources=[
            VolumetricSource(
                block_ids=[source_block_id], group_strength=[source_strength / source_area]
            )
        ]
    )
    forward_solver = SteadyStateSourceSolver(problem=forward_problem)
    forward_solver.Initialize()
    forward_solver.Execute()
    interpolator = FieldFunctionInterpolationVolume()
    interpolator.SetOperationType("sum")
    interpolator.SetLogicalVolume(detector_volume)
    interpolator.AddFieldFunction(forward_problem.GetScalarFluxFieldFunction()[0])
    interpolator.Execute()
    forward_qoi = detector_sigma_a * interpolator.GetValue()
    forward_problem.WriteSurfaceAngularFluxes(
        forward_surface_prefix, interior_surfaces=interior_surfaces
    )

    # Adjoint solve that writes restart data, including the angular flux.
    adjoint_problem = make_adjoint_problem(
        {
            "restart_writes_enabled": True,
            "write_angular_flux_to_restart": True,
            "write_restart_path": restart_stem,
        }
    )
    adjoint_solver = SteadyStateSourceSolver(problem=adjoint_problem)
    adjoint_solver.Initialize()
    adjoint_solver.Execute()
    adjoint_problem.WriteSurfaceAngularFluxes(
        adjoint_surface_prefix, interior_surfaces=interior_surfaces
    )

    # A new adjoint problem reads the restart data when its solver is initialized. Its solution is
    # then in the internal sweep orientation, so writing surface data must be rejected.
    restarted_problem = make_adjoint_problem({"read_restart_path": restart_stem})
    restarted_solver = SteadyStateSourceSolver(problem=restarted_problem)
    restarted_solver.Initialize()
    internal_orientation_rejected = False
    try:
        restarted_problem.WriteSurfaceAngularFluxes(
            rejected_surface_prefix, interior_surfaces=interior_surfaces
        )
    except ValueError as error:
        internal_orientation_rejected = "physical orientation" in str(error)
    if not internal_orientation_rejected:
        raise RuntimeError("Adjoint surface export from restart data was not rejected")
    if os.path.exists(f"{rejected_surface_prefix}{rank}.h5"):
        raise RuntimeError("Rejected adjoint surface export created an output file")

    # Executing the adjoint solver after the restart read leaves the physical orientation.
    restarted_solver.Execute()
    restarted_problem.WriteSurfaceAngularFluxes(
        restarted_surface_prefix, interior_surfaces=interior_surfaces
    )

    def surface_data(name, prefix):
        return {
            "name": name,
            "file_prefixes": {"surface_angular_fluxes": prefix},
            "surfaces": interior_tags,
        }

    evaluator = ResponseEvaluator(problem=forward_problem)
    evaluator.SetOptions(
        data=[
            surface_data("forward", forward_surface_prefix),
            surface_data("adjoint", adjoint_surface_prefix),
            surface_data("restarted_adjoint", restarted_surface_prefix),
        ]
    )
    (surface,) = evaluator.EvaluateSurfaceResponse("forward", "adjoint").values()
    (restarted_surface,) = evaluator.EvaluateSurfaceResponse(
        "forward", "restarted_adjoint"
    ).values()
    surface_qoi = surface["total"]
    restarted_surface_qoi = restarted_surface["total"]

    if rank == 0:
        print(f"Forward QoI={forward_qoi:.8e}")
        print(f"Surface QoI={surface_qoi:.8e}")
        print(f"Restarted QoI={restarted_surface_qoi:.8e}")
        # Relative errors of the surface identity: each surface QoI must equal the forward QoI.
        print(f"Surface identity error={abs(surface_qoi / forward_qoi - 1.0):.3e}")
        print(f"Restarted identity error={abs(restarted_surface_qoi / forward_qoi - 1.0):.3e}")

    barrier()
    os.remove(f"{restart_stem}{rank}.restart.h5")
    for file_prefix in (forward_surface_prefix, adjoint_surface_prefix, restarted_surface_prefix):
        os.remove(f"{file_prefix}{rank}.h5")
