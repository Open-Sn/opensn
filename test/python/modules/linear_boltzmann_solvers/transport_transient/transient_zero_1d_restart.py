#!/usr/bin/env python3

import os
import shutil
import sys

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

if "opensn_console" not in globals():
    sys.path.append("/Users/dhawkins/opensn_dev/opensn")
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver


def write_xs_file(xs_path):
    if rank == 0:
        with open(xs_path, "w", encoding="ascii") as xs_file:
            xs_file.write(
                """NUM_GROUPS 2
NUM_MOMENTS 1

INV_VELOCITY_BEGIN
0 1.0
1 2.0
INV_VELOCITY_END

SIGMA_T_BEGIN
0 1.0
1 1.2
SIGMA_T_END

SIGMA_A_BEGIN
0 0.45
1 0.55
SIGMA_A_END

TRANSFER_MOMENTS_BEGIN
M_GFROM_GTO_VAL 0 0 0 0.35
M_GFROM_GTO_VAL 0 0 1 0.10
M_GFROM_GTO_VAL 0 1 0 0.20
M_GFROM_GTO_VAL 0 1 1 0.35
TRANSFER_MOMENTS_END
"""
            )
    comm.Barrier()


def make_problem(xs_path, restart_read="", restart_write=""):
    dx = 1.0 / 20
    nodes = [i * dx for i in range(21)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_path)

    source = VolumetricSource(
        block_ids=[0], group_strength=[1.5, 0.5], start_time=0.0, end_time=0.30
    )
    quad = GLProductQuadrature1DSlab(n_polar=2, scattering_order=0)

    options = {
        "save_angular_flux": True,
        "verbose_inner_iterations": True,
        "verbose_outer_iterations": False,
    }
    if restart_read:
        options["read_restart_path"] = restart_read
    if restart_write:
        options["restart_writes_enabled"] = True
        options["write_restart_path"] = restart_write

    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        time_dependent=True,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": quad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options=options,
    )


def max_abs_diff(v1, v2):
    diff = 0.0
    for a, b in zip(v1, v2):
        diff = max(diff, abs(a - b))
    return diff


if __name__ == "__main__":
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    restart_dir = "transport_transient_1d_restart"
    restart_stem = os.path.join(restart_dir, "transient")
    xs_path = os.path.join(restart_dir, "two_group_td.xs")
    if rank == 0 and os.path.isdir(restart_dir):
        shutil.rmtree(restart_dir)
    if rank == 0:
        os.makedirs(restart_dir, exist_ok=True)
    comm.Barrier()
    write_xs_file(xs_path)

    dt = 0.05
    split_time = 2.0 * dt
    stop_time = 4.0 * dt

    ref_problem = make_problem(xs_path)
    ref_solver = TransientSolver(
        problem=ref_problem,
        dt=dt,
        theta=1.0,
        stop_time=stop_time,
        initial_state="zero",
        verbose=True,
    )
    ref_solver.Initialize()
    ref_solver.Execute()
    ref_phi = list(ref_problem.GetPhiNewLocal())

    write_problem = make_problem(xs_path, restart_write=restart_stem)
    write_solver = TransientSolver(
        problem=write_problem,
        dt=dt,
        theta=1.0,
        stop_time=split_time,
        initial_state="zero",
        verbose=True,
    )
    write_solver.Initialize()
    write_solver.Execute()

    restart_problem = make_problem(xs_path, restart_read=restart_stem)
    restart_solver = TransientSolver(
        problem=restart_problem,
        dt=dt,
        theta=1.0,
        stop_time=stop_time,
        initial_state="existing",
        verbose=True,
    )
    restart_solver.Initialize()
    restart_solver.Execute()
    restart_phi = list(restart_problem.GetPhiNewLocal())

    phi_diff = max_abs_diff(ref_phi, restart_phi)
    pass_flag = int(
        abs(restart_problem.GetTime() - stop_time) < 1.0e-12 and phi_diff < 1.0e-10
    )

    if rank == 0:
        print(f"TRANSIENT_RESTART_MAX_ABS_DIFF {phi_diff:.16e}")
        print(f"TRANSIENT_RESTART_PASS {pass_flag}")

    comm.Barrier()
    if rank == 0 and os.path.isdir(restart_dir):
        shutil.rmtree(restart_dir)
