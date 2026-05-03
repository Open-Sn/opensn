#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Time-dependent boundary-condition regression tests.

The cases use small absorber problems and compare equivalent boundary formulations where possible.
Zero-inflow and validation cases provide analytic checks for the inactive-boundary behavior and
input.
"""

import math
import os
import sys

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.size
rank = comm.rank

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.math import AngularFluxFunction, AngularFluxTimeFunction
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver


DT = 0.1
TOL = 1.0e-10
XS2_PATH = "td_boundary_2g_tmp.xs"


def ensure_multigroup_xs_file():
    if rank == 0:
        with open(XS2_PATH, "w", encoding="ascii") as xs_file:
            xs_file.write(
                """NUM_GROUPS 2
NUM_MOMENTS 1

INV_VELOCITY_BEGIN
0 1.0
1 1.0
INV_VELOCITY_END

SIGMA_T_BEGIN
0 1.0
1 1.15
SIGMA_T_END

SIGMA_A_BEGIN
0 1.0
1 1.15
SIGMA_A_END

TRANSFER_MOMENTS_BEGIN
TRANSFER_MOMENTS_END
"""
            )
    comm.Barrier()


def make_1d_problem(boundary_conditions, num_groups=1, n_polar=4):
    nodes = [i / 10 for i in range(11)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    if num_groups == 1:
        xs.CreateSimpleOneGroup(1.0, 0.0, 1.0)
    else:
        ensure_multigroup_xs_file()
        xs.LoadFromOpenSn(XS2_PATH)

    pquad = GLProductQuadrature1DSlab(n_polar=n_polar, scattering_order=0)
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        time_dependent=True,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=boundary_conditions,
        options={
            "save_angular_flux": True,
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )


def run_steps(problem, steps, dt=DT):
    solver = TransientSolver(problem=problem, dt=dt, theta=1.0, initial_state="zero", verbose=False)
    solver.Initialize()
    solver.SetTimeStep(dt)
    for _ in range(steps):
        solver.Advance()
    return list(problem.GetPhiNewLocal())


def max_abs(v):
    return max([abs(x) for x in v] + [0.0])


def max_abs_diff(v1, v2):
    return max([abs(a - b) for a, b in zip(v1, v2)] + [0.0])


def global_max(value):
    return comm.allreduce(value, op=MPI.MAX)


def compare_problem_pair(problem_a, problem_b, steps):
    phi_a = run_steps(problem_a, steps)
    phi_b = run_steps(problem_b, steps)
    return global_max(max_abs_diff(phi_a, phi_b))


def run_zero_check(problem, steps):
    return global_max(max_abs(run_steps(problem, steps)))


def iso_bc(name, strength, start=-math.inf, end=math.inf):
    return {
        "name": name,
        "type": "isotropic",
        "group_strength": strength,
        "start_time": start,
        "end_time": end,
    }


def arb_time_bc(name, func):
    return {"name": name, "type": "arbitrary", "time_function": AngularFluxTimeFunction(func)}


def arb_const_bc(name, func):
    return {"name": name, "type": "arbitrary", "function": AngularFluxFunction(func)}


def vacuum_bc(name):
    return {"name": name, "type": "vacuum"}


def reflecting_bc(name):
    return {"name": name, "type": "reflecting"}


def report(name, pass_flag, detail):
    if rank == 0:
        print(f"{name}_DETAIL {detail:.16e}")
        print(f"{name}_PASS {int(pass_flag)}")


def run_case_01():
    strength = 2.0
    diff = compare_problem_pair(
        make_1d_problem([iso_bc("zmin", [strength], 0.0, 0.0), vacuum_bc("zmax")]),
        make_1d_problem(
            [arb_time_bc("zmin", lambda g, n, t: strength if t <= 0.0 else 0.0), vacuum_bc("zmax")]
        ),
        2,
    )
    report("TD_BOUNDARY_CASE01_ISO_PULSE_EQUIV", diff < TOL, diff)


def run_case_02():
    strength = 1.4
    diff = compare_problem_pair(
        make_1d_problem([iso_bc("zmin", [strength], DT, DT), vacuum_bc("zmax")]),
        make_1d_problem(
            [
                arb_time_bc("zmin", lambda g, n, t: strength if abs(t - DT) < 1.0e-12 else 0.0),
                vacuum_bc("zmax"),
            ]
        ),
        3,
    )
    report("TD_BOUNDARY_CASE02_DELAYED_ISO_EQUIV", diff < TOL, diff)


def run_case_03():
    strength = 1.7
    diff = compare_problem_pair(
        make_1d_problem([iso_bc("zmin", [strength], 0.0, DT), vacuum_bc("zmax")]),
        make_1d_problem(
            [
                arb_time_bc("zmin", lambda g, n, t: strength if 0.0 <= t <= DT else 0.0),
                vacuum_bc("zmax"),
            ]
        ),
        3,
    )
    report("TD_BOUNDARY_CASE03_ISO_WINDOW_EQUIV", diff < TOL, diff)


def run_case_04():
    residual = run_zero_check(
        make_1d_problem([iso_bc("zmin", [3.0], 1.0, 2.0), vacuum_bc("zmax")]), 3
    )
    report("TD_BOUNDARY_CASE04_INACTIVE_ISO_ZERO", residual < TOL, residual)


def run_case_05():
    residual = run_zero_check(
        make_1d_problem(
            [arb_time_bc("zmin", lambda g, n, t: 4.0 if t >= 1.0 else 0.0), vacuum_bc("zmax")]
        ),
        3,
    )
    report("TD_BOUNDARY_CASE05_INACTIVE_ARBITRARY_ZERO", residual < TOL, residual)


def run_case_06():
    observed = set()

    def record_time(group, direction, time):
        observed.add(round(time, 12))
        return 0.25 + time

    problem = make_1d_problem([arb_time_bc("zmin", record_time), vacuum_bc("zmax")])
    run_steps(problem, 3)

    local_mask = 0
    for i, t in enumerate([0.0, DT, 2.0 * DT]):
        if round(t, 12) in observed:
            local_mask |= 1 << i
    global_mask = comm.allreduce(local_mask, op=MPI.BOR)
    missing = 0 if global_mask == 0b111 else 1
    report("TD_BOUNDARY_CASE06_ARBITRARY_CALLED_EACH_STEP", missing == 0, float(missing))


def run_case_07():
    strengths = [1.0, 0.35]

    def strength(group, direction, time):
        return strengths[group] if time <= DT else 0.0

    diff = compare_problem_pair(
        make_1d_problem([iso_bc("zmin", strengths, 0.0, DT), vacuum_bc("zmax")], num_groups=2),
        make_1d_problem([arb_time_bc("zmin", strength), vacuum_bc("zmax")], num_groups=2),
        3,
    )
    report("TD_BOUNDARY_CASE07_MULTIGROUP_EQUIV", diff < TOL, diff)


def run_case_08():
    def left(group, direction, time):
        return 1.2 if time <= DT else 0.0

    def right_const(group, direction):
        return 0.4

    def right_time(group, direction, time):
        return 0.4

    diff = compare_problem_pair(
        make_1d_problem([arb_time_bc("zmin", left), arb_const_bc("zmax", right_const)]),
        make_1d_problem([arb_time_bc("zmin", left), arb_time_bc("zmax", right_time)]),
        3,
    )
    report("TD_BOUNDARY_CASE08_MULTI_ARBITRARY_CONSTANT_EQUIV", diff < TOL, diff)


def run_case_09():
    def left(group, direction, time):
        return 1.0 if time <= 0.0 else 0.0

    def right(group, direction, time):
        return 0.6 if abs(time - DT) < 1.0e-12 else 0.0

    def left_ref(group, direction, time):
        return 1.0 if time <= 0.0 else 0.0

    def right_ref(group, direction, time):
        return 0.6 if abs(time - DT) < 1.0e-12 else 0.0

    diff = compare_problem_pair(
        make_1d_problem([arb_time_bc("zmin", left), arb_time_bc("zmax", right)]),
        make_1d_problem([arb_time_bc("zmin", left_ref), arb_time_bc("zmax", right_ref)]),
        3,
    )
    report("TD_BOUNDARY_CASE09_MULTI_ARBITRARY_WINDOWS", diff < TOL, diff)


def run_case_10():
    iso_strength = 0.8
    arb_strength = 0.45
    combo_bcs = [
        reflecting_bc("xmin"),
        reflecting_bc("xmax"),
        arb_time_bc("ymin", lambda g, n, t: arb_strength if t <= DT else 0.0),
        iso_bc("ymax", [iso_strength], 0.0, DT),
        vacuum_bc("zmin"),
        vacuum_bc("zmax"),
    ]
    arbitrary_ok = (
        combo_bcs[2]["time_function"](0, 0, 0.0) == arb_strength
        and combo_bcs[2]["time_function"](0, 0, 2.0 * DT) == 0.0
    )
    isotropic_ok = (
        combo_bcs[3]["group_strength"][0] == iso_strength
        and combo_bcs[3]["start_time"] == 0.0
        and combo_bcs[3]["end_time"] == DT
    )
    reflecting_ok = combo_bcs[0]["type"] == "reflecting" and combo_bcs[1]["type"] == "reflecting"
    pass_flag = arbitrary_ok and isotropic_ok and reflecting_ok
    detail = 0.0 if pass_flag else 1.0
    report("TD_BOUNDARY_CASE10_REFLECTING_ARBITRARY_ISO", pass_flag, detail)


def run_case_11():
    failed_as_expected = False
    try:
        make_1d_problem([iso_bc("zmin", [1.0], DT, 0.0), vacuum_bc("zmax")])
    except Exception:
        failed_as_expected = True
    detail = 0.0 if failed_as_expected else 1.0
    report("TD_BOUNDARY_CASE11_INVALID_ISO_BOUNDS", failed_as_expected, detail)


def run_case_12():
    failed_as_expected = False
    try:
        make_1d_problem(
            [
                {
                    "name": "zmin",
                    "type": "arbitrary",
                    "time_function": AngularFluxTimeFunction(lambda g, n, t: 1.0),
                    "start_time": 0.0,
                },
                vacuum_bc("zmax"),
            ]
        )
    except Exception:
        failed_as_expected = True
    detail = 0.0 if failed_as_expected else 1.0
    report("TD_BOUNDARY_CASE12_ARBITRARY_REJECTS_BOUNDS", failed_as_expected, detail)


if __name__ == "__main__":
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    run_case_01()
    run_case_02()
    run_case_03()
    run_case_04()
    run_case_05()
    run_case_06()
    run_case_07()
    run_case_08()
    run_case_09()
    run_case_10()
    run_case_11()
    run_case_12()

    comm.Barrier()
    if rank == 0 and os.path.exists(XS2_PATH):
        os.remove(XS2_PATH)
