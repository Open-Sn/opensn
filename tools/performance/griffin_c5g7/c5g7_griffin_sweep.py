#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""2D C5G7 sweep benchmark matching the Griffin GPU sweep paper.

N. Choi, J. Hanophy, and Y. Wang, Nucl. Sci. Eng. (2026),
https://doi.org/10.1080/00295639.2026.2700698, Sec. IV.A, Figs. 6, 9, and 10.

This is a performance input, not a regression test. The paper setup it reproduces:
101,728 quadrilaterals (see make_c5g7_quad_mesh.py), first-order DFEM (4 DOFs per quad,
the same as PWLD on quads), 7 groups, P0 scattering, a single domain, and 24 to 384
directions.

The comparison metric is OpenSn's per-sweep log line
    WGS groups [0-6] avg_sweep_time = ... s, sweep_time_per_unknown = ... ns
whose normalization is sweep time / (nodes * angles * groups). The paper defines grind time
the same way (sweep time / (elements * 4 DOFs * 7 groups * directions)).

Parameters (pass with -p name=value):
  num_dirs          swept directions: 24, 48, 96, 192, or 384 (default 384)
  n_polar           GLC polar angles over the full sphere; 2D uses half (default 4)
  mode              "timing" (fixed sweep count) or "keff" (converged k with CMFD)
  use_gpus          default True
  sweep_type        "AAH" or "CBC" (default "AAH")
  angle_agg         "polar", "single", or "azimuthal" (default "polar")
  sweeps_per_outer  timing mode: source-iteration sweeps per power iteration (default 10)
  num_outers        timing mode: power iterations (default 3; the first is warm-up)
  mesh_file         default: c5g7_griffin.msh next to this script
  vacuum_boundaries replace the reflecting xmin/ymin boundaries with vacuum (default False).
                    Diagnostic only: removes the reflecting dependencies between angle sets,
                    so every angle set can sweep at once, without changing the sweep work.
  extra_problem_params  dict merged into the DiscreteOrdinatesProblem arguments
  extra_options     dict merged into the problem's options block, for example
                    {'fused_aao_sweep': True, 'fused_aao_lanes': 8}

Timing mode does not save the angular flux. The paper's GPU sweeper stored the full angular
flux on the device, while OpenSn's GPU path would copy it back to the host after every sweep.
The paper reports that wavefront storage (no full angular flux) changes the P0 sweep time by
at most 4%.

Example:
  python3 make_c5g7_quad_mesh.py c5g7_griffin.msh
  opensn -i c5g7_griffin_sweep.py -p num_dirs=384 -p sweep_type='"AAH"'
"""

import os

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration


def get_option(name, default):
    return globals()[name] if name in globals() else default


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, "../../.."))
    xs_dir = os.path.join(repo_root, "test/assets/xs/c5g7")

    num_dirs = int(get_option("num_dirs", 384))
    n_polar = int(get_option("n_polar", 4))
    mode = get_option("mode", "timing")
    use_gpus = bool(get_option("use_gpus", True))
    sweep_type = get_option("sweep_type", "AAH")
    angle_agg = get_option("angle_agg", "polar")
    sweeps_per_outer = int(get_option("sweeps_per_outer", 10))
    num_outers = int(get_option("num_outers", 3))
    mesh_file = get_option("mesh_file", os.path.join(script_dir, "c5g7_griffin.msh"))
    vacuum_boundaries = bool(get_option("vacuum_boundaries", False))
    extra_problem_params = dict(get_option("extra_problem_params", {}))
    extra_options = dict(get_option("extra_options", {}))

    if mode not in ("timing", "keff"):
        raise ValueError("mode must be 'timing' or 'keff'")
    if n_polar % 2 != 0 or num_dirs % (n_polar // 2) != 0:
        raise ValueError(f"num_dirs={num_dirs} is not divisible by n_polar/2={n_polar // 2}")
    n_azimuthal = num_dirs // (n_polar // 2)
    if n_azimuthal % 4 != 0:
        raise ValueError(f"n_azimuthal={n_azimuthal} must be a multiple of 4")
    if not os.path.isfile(mesh_file):
        raise FileNotFoundError(
            f"{mesh_file} not found; generate it with make_c5g7_quad_mesh.py")

    grid = FromFileMeshGenerator(filename=mesh_file).Execute()
    grid.SetOrthogonalBoundaries()

    xs_files = ["XS_water.xs", "XS_UO2.xs", "XS_7pMOX.xs", "XS_guide_tube.xs",
                "XS_4_3pMOX.xs", "XS_8_7pMOX.xs", "XS_fission_chamber.xs"]
    xs_map = []
    for block_id, xs_file in enumerate(xs_files):
        xs = MultiGroupXS()
        xs.LoadFromOpenSn(os.path.join(xs_dir, xs_file))
        xs_map.append({"block_ids": [block_id], "xs": xs})
    num_groups = xs_map[0]["xs"].num_groups

    pquad = GLCProductQuadrature2DXY(n_polar=n_polar, n_azimuthal=n_azimuthal,
                                     scattering_order=0)

    if mode == "timing":
        # Tolerances that cannot be met force a fixed number of sweeps per outer.
        groupset_tol = {"inner_linear_method": "classic_richardson",
                        "l_max_its": sweeps_per_outer, "l_abs_tol": 1.0e-18}
    else:
        groupset_tol = {"inner_linear_method": "petsc_gmres",
                        "l_max_its": 5, "l_abs_tol": 1.0e-10}

    problem_params = dict(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[{"groups_from_to": (0, num_groups - 1),
                    "angular_quadrature": pquad,
                    "angle_aggregation_type": angle_agg,
                    **groupset_tol}],
        xs_map=xs_map,
        boundary_conditions=[] if vacuum_boundaries else [
            {"name": "xmin", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"}],
        options={"verbose_inner_iterations": True,
                 "verbose_outer_iterations": True,
                 "save_angular_flux": False,
                 "power_default_kappa": 1.0,
                 **extra_options},
        sweep_type=sweep_type,
        use_gpus=use_gpus,
    )
    problem_params.update(extra_problem_params)

    if rank == 0:
        print(f"GRIFFIN_C5G7 config: mode={mode} num_dirs={num_dirs} n_polar={n_polar} "
              f"n_azimuthal={n_azimuthal} sweep_type={sweep_type} angle_agg={angle_agg} "
              f"use_gpus={use_gpus} vacuum_boundaries={vacuum_boundaries} ranks={size} "
              f"sweeps_per_outer={sweeps_per_outer} "
              f"num_outers={num_outers} extra={extra_problem_params} "
              f"extra_options={extra_options}")

    phys = DiscreteOrdinatesProblem(**problem_params)

    if mode == "timing":
        k_solver = PowerIterationKEigenSolver(problem=phys, k_tol=1.0e-18, max_iters=num_outers)
    else:
        # Same CMFD settings as test/python/.../transport_keigen/c5g7/c5g7.py.
        cmfd = CMFDAcceleration(problem=phys, coarse_mesh="local_aggregation",
                                aggregation_size=16, update_wgs_max_its=4,
                                update_wgs_abs_tol=1.0e-4, l_abs_tol=1.0e-10, max_iters=100,
                                pi_max_its=5, pi_k_tol=1.0e-10)
        k_solver = PowerIterationKEigenSolver(problem=phys, acceleration=cmfd, k_tol=1.0e-8)

    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()
    if rank == 0:
        print(f"GRIFFIN_C5G7 k-eigenvalue: {k:.8f}")
