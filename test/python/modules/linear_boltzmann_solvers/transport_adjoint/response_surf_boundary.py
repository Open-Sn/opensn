#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Boundary surface responses and the discrete adjoint identity, in 2D (one group) or 3D (three
groups with upscatter, in two groupsets); option ``case`` = "2d" or "3d".

For the upwind DFEM discretization, the forward and adjoint solutions satisfy

    <q, psi^dagger> - <q^dagger, psi> = sum over boundary faces of
        sum_d w_d mu_d M (psi_up psi^dagger_up),

where on a boundary face, for outflow directions (mu_d > 0), psi_up is the forward cell value and
psi^dagger_up the adjoint boundary-condition value, upwind for the adjoint because it travels along
-Omega_d, and for inflow directions psi_up is the forward boundary-condition value and
psi^dagger_up the adjoint cell value. Interior faces cancel in pairs, and the terms of a reflecting
boundary cancel when both problems reflect. The forward problem of A and B has a volumetric source q
near the origin, vacuum boundaries, and a reflecting ymin boundary. A and B are checked against
<q, phi^dagger> from EvaluateResponse:

A. Adjoint driven by a volumetric detector source with vacuum boundaries (ymin reflecting): the
   interior plane x = x_p, which separates source and detector, carries the response, and the
   vacuum boundary responses are exactly zero (the adjoint boundary values are zero).
B. Adjoint driven only by a direction- and group-dependent incoming flux on xmax, a detector outside
   the domain: the xmax boundary response and the interior plane response both equal the response.
   The boundary values are those of the adjoint sweep, so this checks the mapping of adjoint sweep
   directions to physical directions.

C. Forward problem driven only by an isotropic incoming flux on xmin, with the adjoint of A: the
   xmin incoming part (forward boundary value against the adjoint cell value) equals minus the
   detector response computed from the forward scalar flux, the reflecting ymin terms cancel, and
   the plane carries the response.

It also checks that the upwind values of the reflecting boundary are the reflected angular fluxes.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    barrier = MPI.COMM_WORLD.Barrier

    def global_max(value):
        return MPI.COMM_WORLD.allreduce(value, op=MPI.MAX)

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY, GLCProductQuadrature3DXYZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import AngularFluxFunction
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.response import ResponseEvaluator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier

    def global_max(value):
        return MPIAllReduce(value, "max")


def get_option(name, default):
    return globals()[name] if name in globals() else default


if __name__ == "__main__":
    case = get_option("case", "2d")
    three_d = case == "3d"
    num_procs = 4 if three_d else 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Mesh: source block at the origin, detector block in the far corner, plane x = x_p between.
    if three_d:
        length, num_cells, plane = 4.0, 8, 2.5
        source_max, detector_min = 1.0, 3.0
    else:
        length, num_cells, plane = 10.0, 20, 6.0
        source_max, detector_min = 2.0, 8.0
    nodes = [i * length / num_cells for i in range(num_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes] * (3 if three_d else 2)).Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    infz = not three_d
    zlim = {} if infz else {"zmin": 0.0, "zmax": source_max}
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(
            xmin=0.0, xmax=source_max, ymin=0.0, ymax=source_max, infz=infz, **zlim
        ),
        1,
        True,
    )
    zlim = {} if infz else {"zmin": detector_min, "zmax": length}
    detector_volume = RPPLogicalVolume(
        xmin=detector_min, xmax=length, ymin=detector_min, ymax=length, infz=infz, **zlim
    )
    grid.SetBlockIDFromLogicalVolume(detector_volume, 2, True)

    if three_d:
        num_groups = 3
        xs = MultiGroupXS()
        xs.LoadFromOpenSn("../../../../assets/xs/simple_upscatter.xs")
        xs_map = [{"block_ids": [0, 1, 2], "xs": xs}]
        quadrature = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=0)
        group_ranges = [(0, 0), (1, 2)]
        boundaries = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]
    else:
        num_groups = 1
        xs_background = MultiGroupXS()
        xs_background.CreateSimpleOneGroup(0.3, 0.8)
        xs_source = MultiGroupXS()
        xs_source.CreateSimpleOneGroup(1.0, 0.5)
        xs_detector = MultiGroupXS()
        xs_detector.CreateSimpleOneGroup(0.8, 0.0)
        xs_map = [
            {"block_ids": [0], "xs": xs_background},
            {"block_ids": [1], "xs": xs_source},
            {"block_ids": [2], "xs": xs_detector},
        ]
        quadrature = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16, scattering_order=0)
        group_ranges = [(0, 0)]
        boundaries = ["xmin", "xmax", "ymin", "ymax"]
    vacuum_boundaries = [b for b in boundaries if b != "ymin"]

    def make_problem(adjoint, xmax_condition, xmin_condition=None):
        conditions = [{"name": "ymin", "type": "reflecting"}, xmax_condition]
        if xmin_condition is not None:
            conditions.append(xmin_condition)
        conditions += [
            {"name": b, "type": "vacuum"}
            for b in vacuum_boundaries
            if b != "xmax" and not (xmin_condition is not None and b == "xmin")
        ]
        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=num_groups,
            groupsets=[
                {
                    "groups_from_to": groups,
                    "angular_quadrature": quadrature,
                    "inner_linear_method": "petsc_gmres",
                    "l_abs_tol": 1.0e-13,
                    "l_max_its": 1000,
                    "gmres_restart_interval": 100,
                }
                for groups in group_ranges
            ],
            xs_map=xs_map,
            boundary_conditions=conditions,
            options={
                "save_angular_flux": True,
                "adjoint": adjoint,
                "max_ags_iterations": 500,
                "ags_tolerance": 1.0e-13,
            },
        )

    def solve(problem):
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()

    interior = {"plane": ("x", plane)}
    prefix = f"SurfBoundary{case}_"
    vacuum = {"name": "xmax", "type": "vacuum"}

    source_strength = [1.0] * num_groups
    forward = make_problem(False, vacuum)
    forward.SetVolumetricSources(
        volumetric_sources=[VolumetricSource(block_ids=[1], group_strength=source_strength)]
    )
    solve(forward)
    forward.WriteSurfaceAngularFluxes(
        prefix + "F_", boundary_surfaces=boundaries, interior_surfaces=interior
    )

    # A. Volumetric detector source, vacuum adjoint boundaries.
    adjoint_a = make_problem(True, vacuum)
    adjoint_a.SetVolumetricSources(
        volumetric_sources=[
            VolumetricSource(logical_volume=detector_volume, group_strength=[0.8] * num_groups)
        ]
    )
    solve(adjoint_a)
    adjoint_a.WriteFluxMoments(prefix + "A_phi_")
    adjoint_a.WriteSurfaceAngularFluxes(
        prefix + "A_", boundary_surfaces=boundaries, interior_surfaces=interior
    )

    # B. Incoming adjoint flux on xmax that depends on the direction and the group.
    def incoming(group, direction):
        return 1.0 + 0.1 * direction + 0.3 * group

    adjoint_b = make_problem(
        True, {"name": "xmax", "type": "arbitrary", "function": AngularFluxFunction(incoming)}
    )
    solve(adjoint_b)
    adjoint_b.WriteFluxMoments(prefix + "B_phi_")
    adjoint_b.WriteSurfaceAngularFluxes(
        prefix + "B_", boundary_surfaces=boundaries, interior_surfaces=interior
    )

    # C. Forward problem driven only by an isotropic incoming flux on xmin; its detector response,
    # sum over groups of 0.8 times the integral of the scalar flux over the detector, is computed
    # from the forward flux independently of the surface data.
    forward_c = make_problem(
        False, vacuum, {"name": "xmin", "type": "isotropic", "group_strength": [1.0] * num_groups}
    )
    solve(forward_c)
    qoi_c = 0.0
    for field_function in forward_c.GetScalarFluxFieldFunction():
        interpolator = FieldFunctionInterpolationVolume()
        interpolator.SetOperationType("sum")
        interpolator.SetLogicalVolume(detector_volume)
        interpolator.AddFieldFunction(field_function)
        interpolator.Execute()
        qoi_c += 0.8 * interpolator.GetValue()
    forward_c.WriteSurfaceAngularFluxes(
        prefix + "C_", boundary_surfaces=boundaries, interior_surfaces=interior
    )

    surfaces = boundaries + ["plane"]

    def data(name, flux_prefix, surface_prefix):
        file_prefixes = {"surface_angular_fluxes": surface_prefix}
        if flux_prefix:
            file_prefixes["flux_moments"] = flux_prefix
        return {"name": name, "file_prefixes": file_prefixes, "surfaces": surfaces}

    evaluator = ResponseEvaluator(problem=forward)
    evaluator.SetOptions(
        data=[
            data("F", None, prefix + "F_"),
            data("A", prefix + "A_phi_", prefix + "A_"),
            data("B", prefix + "B_phi_", prefix + "B_"),
            data("C", None, prefix + "C_"),
        ],
        sources={"material": [{"block_id": 1, "strength": source_strength}]},
    )
    qoi_a = evaluator.EvaluateResponse("A")
    qoi_b = evaluator.EvaluateResponse("B")
    plane_a = evaluator.EvaluateSurfaceResponse("F", "A", ["plane"])["plane"]["total"]
    vacuum_a = sum(
        r["total"] for r in evaluator.EvaluateSurfaceResponse("F", "A", vacuum_boundaries).values()
    )
    plane_b = evaluator.EvaluateSurfaceResponse("F", "B", ["plane"])["plane"]["total"]
    xmax_b = evaluator.EvaluateSurfaceResponse("F", "B", ["xmax"])["xmax"]["total"]

    # With no volumetric forward source, the boundary totals sum to minus the detector response:
    # the xmin incoming part carries it (forward boundary value against the adjoint cell value),
    # the vacuum boundaries carry nothing, and the reflecting ymin terms cancel. The plane, which
    # separates xmin from the detector, carries the response.
    responses_c = evaluator.EvaluateSurfaceResponse("C", "A", boundaries + ["plane"])
    xmin_incoming_c = responses_c["xmin"]["incoming"]
    boundary_total_c = sum(responses_c[b]["total"] for b in boundaries)
    ymin_total_c = responses_c["ymin"]["total"]
    plane_c = responses_c["plane"]["total"]

    # A surface list without surface files is rejected.
    stray_surfaces_rejected = 0
    try:
        evaluator.SetOptions(
            data=[{"name": "stray", "file_prefixes": {}, "surfaces": ["xmax"]}]
        )
    except ValueError as error:
        stray_surfaces_rejected = int("surface_angular_fluxes" in str(error))

    # On the reflecting ymin boundary (normal -y), the upwind value of a direction entering the
    # domain is the angular flux of its reflection (omega_y -> -omega_y) at the same node.
    reflection_mismatch = 0.0
    for problem_prefix, adjoint in ((prefix + "F_", False), (prefix + "B_", True)):
        for entry in forward.ReadSurfaceAngularFluxes(problem_prefix, ["ymin"]):
            values = entry["data"]
            mapping = entry["mapping"]
            # Ranks without ymin faces read empty arrays.
            if not values["psi"]:
                continue
            num_dirs = entry["num_directions"]
            groups = entry["num_groups"]
            omega = values["omega"]
            reflected = []
            for d in range(num_dirs):
                target = (omega[3 * d], -omega[3 * d + 1], omega[3 * d + 2])
                reflected.append(
                    min(
                        range(num_dirs),
                        key=lambda r: sum((omega[3 * r + c] - target[c]) ** 2 for c in range(3)),
                    )
                )
            # The face of each node, for its normal.
            node_faces = [
                f for f, n in enumerate(mapping["num_face_nodes"]) for _ in range(n)
            ]
            for node, face in enumerate(node_faces):
                normal = mapping["normals"][3 * face : 3 * face + 3]
                for d in range(num_dirs):
                    mu = sum(omega[3 * d + c] * normal[c] for c in range(3))
                    if not (mu > 0.0 if adjoint else mu < 0.0):
                        continue
                    for g in range(groups):
                        up = values["upwind_psi"][(node * num_dirs + d) * groups + g]
                        own = values["psi"][(node * num_dirs + reflected[d]) * groups + g]
                        reflection_mismatch = max(reflection_mismatch, abs(up - own))
    reflection_mismatch = global_max(reflection_mismatch)

    if rank == 0:
        print(f"A response={qoi_a:.8e}")
        print(f"B response={qoi_b:.8e}")
        print(f"A plane identity error={abs(plane_a / qoi_a - 1.0):.3e}")
        print(f"A vacuum boundary fraction={abs(vacuum_a / qoi_a):.3e}")
        print(f"B plane identity error={abs(plane_b / qoi_b - 1.0):.3e}")
        print(f"B boundary identity error={abs(xmax_b / qoi_b - 1.0):.3e}")
        print(f"C response={qoi_c:.8e}")
        print(f"C xmin incoming identity error={abs(xmin_incoming_c / qoi_c + 1.0):.3e}")
        print(f"C boundary total identity error={abs(boundary_total_c / qoi_c + 1.0):.3e}")
        print(f"C reflecting total fraction={abs(ymin_total_c / qoi_c):.3e}")
        print(f"C plane identity error={abs(plane_c / qoi_c - 1.0):.3e}")
        print(f"Reflected upwind mismatch={reflection_mismatch:.3e}")
        print(f"Stray surfaces rejected={stray_surfaces_rejected}")

    barrier()
    if rank == 0:
        for name in os.listdir("."):
            if name.startswith(prefix):
                os.remove(name)
    barrier()
