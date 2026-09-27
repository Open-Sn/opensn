#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_integrated_quantities.py: RZ integrated quantities in physical units.

A cylinder of radius R = 3 and height H = 2 with a uniform source Q = 1, sigma_t = 1, c = 0.6,
reflecting rmin/zmin/zmax, and isotropic inflow on rmax equal to the infinite-medium flux has
the exact uniform solution phi = Q / sigma_a = 2.5. Every integrated quantity must use the
physical volume 2 pi (R^2 / 2) H = 18 pi:

- balance-table production and absorption rates equal Q * V and sigma_a * phi * V,
- volume field-function integral of phi and VolumePostprocessor integral equal phi * V,
- a response evaluation with a unit material source and phi as the buffer equals phi * V,
- a power field normalized to a specified total integrates to that total, and
- exported finite-element surface weights on rmax integrate to its physical area 2 pi R H.

A separate problem with a point source of strength S = 3 (a ring in RZ) on a mesh vertex shared
by four cells, and vacuum outer
boundaries must report a production rate of exactly S, and ComputeLeakage must accept a single
RZ boundary name and match the balance-table outflow.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    global_sum = MPI.COMM_WORLD.allreduce
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.post import VolumePostprocessor
    from pyopensn.response import ResponseEvaluator
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.source import PointSource, VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier
    global_sum = MPIAllReduce


def make_problem(grid, xs, bcs, **sources):
    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DRZ(
                    n_polar=4, n_azimuthal=8, scattering_order=0
                ),
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=bcs,
        options={
            "save_angular_flux": True,
            "use_precursors": False,
            "verbose_inner_iterations": False,
        },
        **sources,
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return problem, solver


def field_function_volume_integral(field_function, op):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType(op)
    ffi.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    ffi.AddFieldFunction(field_function)
    ffi.Execute()
    return ffi.GetValue()


def volume_integral(problem, op):
    return field_function_volume_integral(problem.GetScalarFluxFieldFunction()[0], op)


if __name__ == "__main__":
    R, H = 3.0, 2.0
    Q, sigma_t, c = 1.0, 1.0, 0.6
    sigma_a = sigma_t * (1.0 - c)
    phi_exact = Q / sigma_a
    volume = 2.0 * math.pi * (R * R / 2.0) * H

    grid = OrthogonalMeshGenerator(
        node_sets=[[R * i / 6 for i in range(7)], [H * i / 4 for i in range(5)]],
        coord_sys="cylindrical",
    ).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, c)

    bcs = [
        {"name": "rmin", "type": "reflecting"},
        {"name": "zmin", "type": "reflecting"},
        {"name": "zmax", "type": "reflecting"},
        {"name": "rmax", "type": "isotropic", "group_strength": [phi_exact]},
    ]
    problem, solver = make_problem(
        grid, xs, bcs, volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[Q])]
    )
    table = solver.ComputeBalanceTable()

    postprocessor = VolumePostprocessor(problem=problem, value_type="integral", group=0)
    postprocessor.Execute()

    prefix = "RZIntegratedPhi_p"
    problem.WriteFluxMoments(prefix)
    evaluator = ResponseEvaluator(problem=problem)
    evaluator.SetOptions(
        buffers=[{"name": "phi", "file_prefixes": {"flux_moments": prefix}}],
        sources={"material": [{"block_id": 0, "strength": [1.0]}]},
    )
    evaluator_integral = evaluator.EvaluateResponse("phi")
    ff_integral = volume_integral(problem, "sum")
    ff_average = volume_integral(problem, "avg")
    vp_integral = float(postprocessor.GetValue()[0][0])

    surface_prefix = "RZIntegratedSurface_p"
    problem.WriteSurfaceAngularFluxes(surface_prefix, boundary_surfaces=["rmax"])
    (surface,) = problem.ReadSurfaceAngularFluxes(surface_prefix, ["rmax"])
    num_surface_nodes = sum(surface["mapping"]["num_face_nodes"])
    local_fe_area = 0.0
    if num_surface_nodes > 0:
        num_directions = len(surface["data"]["fe_shape"]) // num_surface_nodes
        local_fe_area = sum(surface["data"]["fe_shape"]) / num_directions
    fe_surface_area = global_sum(local_fe_area)
    mass_surface_area = global_sum(sum(surface["data"]["M_ij"], 0.0))
    lateral_area = 2.0 * math.pi * R * H

    fissile_xs = MultiGroupXS()
    fissile_xs.LoadFromOpenSn("../../../../assets/xs/simple_fissile_1g.xs")
    problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": fissile_xs}])
    power_target = 7.0
    normalized_power = problem.CreateFieldFunction(
        "normalized_power", "power", power_normalization_target=power_target
    )
    normalized_power_integral = field_function_volume_integral(normalized_power, "sum")

    # Point source (ring) with vacuum outer boundaries.
    S = 3.0
    point_bcs = [{"name": "rmin", "type": "reflecting"}] + [
        {"name": n, "type": "vacuum"} for n in ("rmax", "zmin", "zmax")
    ]
    point_problem, point_solver = make_problem(
        grid, xs, point_bcs, point_sources=[PointSource(location=[1.0, 1.0, 0.0], strength=[S])]
    )
    point_table = point_solver.ComputeBalanceTable()
    leakage = sum(float(point_problem.ComputeLeakage([n])[n][0]) for n in ("rmax", "zmin", "zmax"))
    invalid_rejected = False
    try:
        point_problem.ComputeLeakage(["xmax"])
    except ValueError:
        invalid_rejected = True

    if rank == 0:
        print(f"RZ_PRODUCTION_OVER_QV={table['production_rate'] / (Q * volume):.12e}")
        absorption_ratio = table["absorption_rate"] / (sigma_a * phi_exact * volume)
        print(f"RZ_ABSORPTION_OVER_EXACT={absorption_ratio:.12e}")
        print(f"RZ_FF_SUM_OVER_EXACT={ff_integral / (phi_exact * volume):.12e}")
        print(f"RZ_FF_AVG_OVER_EXACT={ff_average / phi_exact:.12e}")
        print(f"RZ_VOLUME_PP_OVER_EXACT={vp_integral / (phi_exact * volume):.12e}")
        print(f"RZ_EVALUATOR_OVER_EXACT={evaluator_integral / (phi_exact * volume):.12e}")
        print(f"RZ_POWER_NORMALIZATION_OVER_TARGET={normalized_power_integral / power_target:.12e}")
        print(f"RZ_SURFACE_FE_AREA_OVER_EXACT={fe_surface_area / lateral_area:.12e}")
        print(f"RZ_SURFACE_MASS_AREA_OVER_EXACT={mass_surface_area / lateral_area:.12e}")
        print(f"RZ_POINT_SOURCE_PRODUCTION_OVER_S={point_table['production_rate'] / S:.12e}")
        print(f"RZ_LEAKAGE_OVER_OUTFLOW={leakage / point_table['outflow_rate']:.12e}")
        print(f"RZ_INVALID_BOUNDARY_REJECTED={int(invalid_rejected)}")

    barrier()
    os.remove(f"{prefix}{rank}.h5")
    os.remove(f"{surface_prefix}{rank}.h5")
