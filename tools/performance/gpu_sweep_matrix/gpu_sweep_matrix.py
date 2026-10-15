#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fixed-work GPU sweep timing across mesh types, group counts, and angle aggregation.

Performance input, not a regression test. It is meant for A/B comparisons of GPU sweep changes on
problems other than one benchmark: every case runs the same number of source-iteration sweeps of
a steady-state fixed-source problem with vacuum boundaries, so the timing does not depend on
convergence. Compare OpenSn's per-sweep log line

    WGS groups [...] avg_sweep_time = ... s, sweep_time_per_unknown = ... ns

and the printed scalar-flux checksum between builds.

Parameters (pass with -p name=value):
  mesh_kind      "quad2d", "tri2d", "hex3d", "prism3d", or "mixed2d" (default "quad2d")
  num_groups     1, 7, 64, or 168 (default 7)
  n_polar        GLC polar angles over the full sphere (default 4)
  n_azimuthal    GLC azimuthal angles (default 32 in 2D, 16 in 3D; 64 directions either way)
  angle_agg      "polar" or "single" (default "polar")
  sweep_type     "AAH" or "CBC" (default "AAH")
  fused          use the fused AAO GPU scheduler (default True; AAH only)
  sweeps         source-iteration sweeps (default 6)
  work_dir       directory for generated meshes (default: this script's directory)

Meshes (about 65k cells each): quad2d 256x256 quads; tri2d 180x180 squares split into triangles;
hex3d 40^3 hexahedra; prism3d 60x60 triangulated squares extruded into 9 layers; mixed2d the C5G7
coarse mesh (quads and triangles, so mixed cell DOF counts).
"""

import os

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    from pyopensn.mesh import (OrthogonalMeshGenerator, FromFileMeshGenerator,
                               ExtruderMeshGenerator)
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY, GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver


def get_option(name, default):
    return globals()[name] if name in globals() else default


def write_triangle_mesh(path, n, length):
    """Gmsh 2.2 mesh of an n x n grid of squares, each split into two counterclockwise triangles."""
    h = length / n
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write(f"$Nodes\n{(n + 1) ** 2}\n")
        for j in range(n + 1):
            for i in range(n + 1):
                f.write(f"{j * (n + 1) + i + 1} {i * h:.12g} {j * h:.12g} 0\n")
        f.write("$EndNodes\n")
        f.write(f"$Elements\n{2 * n * n}\n")
        e = 0
        for j in range(n):
            for i in range(n):
                a = j * (n + 1) + i + 1
                b, c, d = a + 1, a + n + 2, a + n + 1
                e += 1
                f.write(f"{e} 2 2 0 0 {a} {b} {c}\n")
                e += 1
                f.write(f"{e} 2 2 0 0 {a} {c} {d}\n")
        f.write("$EndElements\n")


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, "../../.."))
    xs_dir = os.path.join(repo_root, "test/assets/xs")

    mesh_kind = get_option("mesh_kind", "quad2d")
    num_groups = int(get_option("num_groups", 7))
    angle_agg = get_option("angle_agg", "polar")
    sweep_type = get_option("sweep_type", "AAH")
    fused = bool(get_option("fused", True)) and sweep_type == "AAH"
    sweeps = int(get_option("sweeps", 6))
    work_dir = get_option("work_dir", script_dir)
    is_3d = mesh_kind in ("hex3d", "prism3d")
    n_polar = int(get_option("n_polar", 4))
    n_azimuthal = int(get_option("n_azimuthal", 16 if is_3d else 32))

    if mesh_kind == "quad2d":
        nodes = [256.0 * i / 256 for i in range(257)]
        grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    elif mesh_kind == "hex3d":
        nodes = [40.0 * i / 40 for i in range(41)]
        grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    elif mesh_kind in ("tri2d", "prism3d"):
        n = 180 if mesh_kind == "tri2d" else 60
        path = os.path.join(work_dir, f"gpu_sweep_matrix_tri{n}.msh")
        if not os.path.isfile(path):
            write_triangle_mesh(path, n, float(n))
        if mesh_kind == "tri2d":
            grid = FromFileMeshGenerator(filename=path).Execute()
        else:
            grid = ExtruderMeshGenerator(
                inputs=[FromFileMeshGenerator(filename=path)],
                layers=[{"z": 9.0, "n": 9}]).Execute()
    elif mesh_kind == "mixed2d":
        grid = FromFileMeshGenerator(
            filename=os.path.join(repo_root, "test/assets/mesh/c5g7/2d_c5g7_coarse.msh")).Execute()
    else:
        raise ValueError(f"unknown mesh_kind {mesh_kind}")
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    if num_groups == 1:
        xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.5)
    elif num_groups == 7:
        xs.LoadFromOpenSn(os.path.join(xs_dir, "c5g7/XS_UO2.xs"))
    elif num_groups == 64:
        xs.LoadFromOpenSn(os.path.join(xs_dir, "diag_XS_64g_1mom_c0.99.xs"))
    elif num_groups == 168:
        xs.LoadFromOpenSn(os.path.join(xs_dir, "xs_168g.xs"))
    else:
        raise ValueError("num_groups must be 1, 7, 64, or 168")
    num_groups = xs.num_groups

    if is_3d:
        pquad = GLCProductQuadrature3DXYZ(n_polar=n_polar, n_azimuthal=n_azimuthal,
                                          scattering_order=0)
    else:
        pquad = GLCProductQuadrature2DXY(n_polar=n_polar, n_azimuthal=n_azimuthal,
                                         scattering_order=0)

    problem_params = dict(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[{"groups_from_to": (0, num_groups - 1),
                    "angular_quadrature": pquad,
                    "angle_aggregation_type": angle_agg,
                    "inner_linear_method": "classic_richardson",
                    "l_max_its": sweeps,
                    "l_abs_tol": 1.0e-18}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0] * num_groups)],
        options={"verbose_inner_iterations": True, "max_ags_iterations": 1},
        sweep_type=sweep_type,
        use_gpus=True,
    )
    if fused:
        problem_params["fused_aao_sweep"] = True

    if rank == 0:
        print(f"GPU_SWEEP_MATRIX config: mesh_kind={mesh_kind} num_groups={num_groups} "
              f"directions={pquad.GetNumAngles() if hasattr(pquad, 'GetNumAngles') else 'n/a'} "
              f"angle_agg={angle_agg} sweep_type={sweep_type} fused={fused} sweeps={sweeps}")

    problem = DiscreteOrdinatesProblem(**problem_params)
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    phi = problem.GetPhiNewLocal()
    if rank == 0:
        print(f"GPU_SWEEP_MATRIX phi_checksum: {sum(phi):.15e}")
