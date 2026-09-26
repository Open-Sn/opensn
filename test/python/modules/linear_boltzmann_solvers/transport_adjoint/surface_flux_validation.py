#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Surface angular-flux export validation.

1. An interior surface whose generated tag (<name>_u or <name>_d) equals a requested boundary
   name is rejected before any output file is created.
2. An interior plane that lies inside the mesh bounds but matches no face, here in the gap
   between two disconnected mesh blocks, is rejected before any output file is created.
3. On a mesh with 1e-8 cm cells translated 1000 cm from the origin, the file reads back and every
   face key on the upper side of an interior surface has a matching key on the lower side.
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
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator, OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier


def make_problem(grid, boundary_names):
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.5)
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=2, n_azimuthal=4, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=[{"name": name, "type": "vacuum"} for name in boundary_names],
        options={"save_angular_flux": True},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return problem


def expect_rejected(problem, prefix, expected_error, **surfaces):
    rejected = False
    try:
        problem.WriteSurfaceAngularFluxes(prefix, **surfaces)
    except ValueError as error:
        rejected = expected_error in str(error)
    if not rejected:
        raise RuntimeError(f"Surface selection {surfaces} was not rejected")
    if os.path.exists(f"{prefix}{rank}.h5"):
        raise RuntimeError("Rejected surface selection created an output file")


def write_two_block_obj(path):
    """Two 2x2 blocks of unit-square quads, x in [0, 1] and [2, 3], y in [0, 1]."""
    vertices = []
    faces = []
    for x0 in (0.0, 2.0):
        base = len(vertices)
        for j in range(3):
            for i in range(3):
                vertices.append((x0 + 0.5 * i, 0.5 * j))
        for j in range(2):
            for i in range(2):
                v = base + 3 * j + i + 1
                faces.append((v, v + 1, v + 4, v + 3))
    with open(path, "w") as f:
        f.write("o TwoBlocks\n")
        for x, y in vertices:
            f.write(f"v {x} {y} 0.0\n")
        for face in faces:
            f.write("f " + " ".join(str(v) for v in face) + "\n")


if __name__ == "__main__":
    num_procs = 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    prefix = "SurfaceValidation_p"
    # Remove output left by an interrupted run so the rejection checks see only new files.
    if os.path.exists(f"{prefix}{rank}.h5"):
        os.remove(f"{prefix}{rank}.h5")
    barrier()

    # 1. A boundary named like a generated interior tag.
    nodes = [0.5 * i for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    grid.SetBoundaryIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.0, xmax=4.0, ymin=3.9, ymax=4.1, infz=True), "plane_u", True
    )
    problem = make_problem(grid, ["xmin", "xmax", "ymin", "plane_u"])
    expect_rejected(
        problem,
        prefix,
        "collides with a requested boundary name",
        boundary_surfaces=["plane_u"],
        interior_surfaces={"plane": ("x", 2.0)},
    )

    # 2. A plane in the gap between two disconnected blocks.
    obj_file = "surface_flux_validation_two_blocks.obj"
    if rank == 0:
        write_two_block_obj(obj_file)
    barrier()
    two_block_grid = FromFileMeshGenerator(filename=obj_file).Execute()
    two_block_grid.SetUniformBlockID(0)
    barrier()
    if rank == 0:
        os.remove(obj_file)
    two_block_problem = make_problem(two_block_grid, [])
    expect_rejected(
        two_block_problem,
        prefix,
        "matches no mesh faces",
        interior_surfaces={"gap": ("x", 1.5)},
    )

    # 3. Face keys on a translated mesh with 1e-8 cm cells.
    scale = 1.0e-8
    origin = 1000.0
    micro_nodes = [origin + scale * i for i in range(9)]
    micro_grid = OrthogonalMeshGenerator(node_sets=[micro_nodes, micro_nodes]).Execute()
    micro_grid.SetOrthogonalBoundaries()
    micro_grid.SetUniformBlockID(0)
    micro_problem = make_problem(micro_grid, ["xmin", "xmax", "ymin", "ymax"])
    micro_problem.WriteSurfaceAngularFluxes(
        prefix, interior_surfaces={"mid": ("x", origin + 4.0 * scale)}
    )
    up, down = micro_problem.ReadSurfaceAngularFluxes(prefix, ["mid_u", "mid_d"])
    up_keys = set(up["mapping"]["cell_map"].keys())
    down_keys = set(down["mapping"]["cell_map"].keys())
    local_faces = len(up["mapping"]["cell_ids"])
    if len(up_keys) != local_faces or len(down_keys) != len(down["mapping"]["cell_ids"]):
        raise RuntimeError("Distinct surface faces share a centroid key")
    # The two-rank partition of this mesh does not cut the x = 4 plane, so both sides of every
    # local face are local and their keys must agree exactly.
    if local_faces == 0 or up_keys != down_keys:
        raise RuntimeError("Upper and lower face keys differ for the same faces")
    barrier()
    os.remove(f"{prefix}{rank}.h5")

    if rank == 0:
        print("SURFACE_VALIDATION_PASSED=1")
