#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Delayed-precursor restarts with two block IDs.

The precursor vectors hold, per node, as many entries as the material with the most precursor
families. The restart files must use this storage stride, not the total family count, which sums
the families over the block IDs of the xs_map. The input variable xs_case selects the xs_map:

  one_material    one 1-group material with one family on block IDs 0 and 1
                  (total 2, stride 1);
  mixed_families  a 1-group material with one family and one with two families, on block IDs
                  0 and 1 in both assignments (total 3, stride 2);
  multigroup      two 2-group materials with one family each, on block IDs 0 and 1
                  (total 2, stride 1);
  header_count    one_material's xs_map (total 2, stride 1) and a 1-group material with two
                  families on block ID 0 beside a non-fissile scatterer on block ID 1 (total 2,
                  stride 2), for the header check described below.

Each xs_map of the first three cases is restarted in two ways, and each restart is compared with
an uninterrupted run that takes the same backward-Euler steps after the same steady state:

  initial condition  a steady-state solve with a uniform source writes a restart; a second problem
                     reads it through read_initial_condition_path and takes one step with the
                     source removed;
  transient restart  after the steady state, the source is removed, one step is taken and a
                     time-dependent restart is written; a second, time-dependent problem reads it
                     through read_restart_path, which also checks the precursor count in the file
                     header, and takes the second step.

The restarts must store and restore the flux and precursor state exactly.

The file header records the total family count, as earlier builds wrote it, and read_restart_path
compares it with the problem's total. Files written with use_precursors=False hold no precursor
data, and they must still read. In header_count, a steady-state restart is written with
use_precursors=False on each xs_map and read through read_restart_path by a problem with
use_precursors=True on the other xs_map, and the flux it reads must equal the written flux. Both
totals are 2. The second xs_map's maximum is also 2, so its header and its check hold 2 whichever
count a build uses: its file stands in for a one_material file written by an earlier build, and
its problem for an earlier build reading a one_material file.

Check:
  PRECURSOR_RESTART_TWO_BLOCKS_DIFF: max over the restarts, phi and psi of
  |split - continuous| / max|continuous phi|; in header_count, max over the two reads of
  |phi read - phi written| / max|phi written|.
"""

import os
import shutil
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op=_mpi_ops[op])

    def MPIBarrier():
        MPI.COMM_WORLD.Barrier()

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, TransientSolver


# For each case: the number of groups and the xs_map assignments, as (block IDs, file) pairs. A
# file None is a non-fissile scatterer.
XS_CASES = {
    "one_material": (1, [[([0, 1], "xs1g_delayed_sub_1p.cxs")]]),
    "mixed_families": (
        1,
        [
            [([0], "xs1g_delayed_sub_1p.cxs"), ([1], "xs1g_delayed_deep_sub_2p_decay_test.cxs")],
            [([0], "xs1g_delayed_deep_sub_2p_decay_test.cxs"), ([1], "xs1g_delayed_sub_1p.cxs")],
        ],
    ),
    "multigroup": (
        2,
        [[([0], "xs2g_delayed_crit_1p.cxs"), ([1], "xs2g_delayed_density_up_1p.cxs")]],
    ),
    "header_count": (
        1,
        [
            [([0, 1], "xs1g_delayed_sub_1p.cxs")],
            [([0], "xs1g_delayed_deep_sub_2p_decay_test.cxs"), ([1], None)],
        ],
    ),
}


def load_xs_map(assignment):
    xs_map = []
    for block_ids, xs_file in assignment:
        xs = MultiGroupXS()
        if xs_file is None:
            xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.5)
        else:
            xs.LoadFromOpenSn(os.path.join("../../../../assets/xs", xs_file))
        xs_map.append({"block_ids": block_ids, "xs": xs})
    return xs_map


def build_problem(grid, num_groups, xs_map, sources, quadrature, options, time_dependent=False):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        time_dependent=time_dependent,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=xs_map,
        volumetric_sources=sources,
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options=options,
    )


def solve_steady_state(phys):
    solver = SteadyStateSourceSolver(problem=phys)
    solver.Initialize()
    solver.Execute()


def run_decay_steps(phys, dt, stop_time, before_execute=None):
    solver = TransientSolver(
        problem=phys, dt=dt, theta=1.0, stop_time=stop_time, initial_state="existing", verbose=False
    )
    solver.Initialize()
    if before_execute is not None:
        before_execute()
    solver.Execute()
    return (
        np.array(phys.GetPhiNewLocal(), copy=True),
        [np.array(psi, copy=True) for psi in phys.GetPsi()],
    )


def global_max_abs(arrays):
    # A NaN counts as an infinite difference.
    local = 0.0
    for a in arrays:
        if a.size > 0:
            value = float(np.max(np.abs(a)))
            local = max(local, np.inf if np.isnan(value) else value)
    return MPIAllReduce(local, "max")


def relative_difference(split, continuous):
    phi_split, psi_split = split
    phi_continuous, psi_continuous = continuous
    if len(psi_split) != len(psi_continuous):
        raise ValueError("angular flux layout differs between the split and continuous runs")
    diff = global_max_abs(
        [phi_split - phi_continuous] + [a - b for a, b in zip(psi_split, psi_continuous)]
    )
    return diff / global_max_abs([phi_continuous])


def restart_differences(setup, ic_base, td_base):
    grid, num_groups, xs_map, sources, quadrature, common_options, dt = setup

    def make_problem(options, problem_sources, time_dependent=False):
        return build_problem(
            grid, num_groups, xs_map, problem_sources, quadrature, options, time_dependent
        )

    def continuous_run(stop_time):
        phys = make_problem(dict(common_options), sources)
        solve_steady_state(phys)
        phys.SetVolumetricSources(clear_volumetric_sources=True)
        phys.SetTimeDependentMode()
        return run_decay_steps(phys, dt, stop_time)

    try:
        # Initial condition: a steady-state restart without angular flux.
        phys_steady = make_problem(
            {
                **common_options,
                "restart_writes_enabled": True,
                "write_angular_flux_to_restart": False,
                "write_delayed_psi_to_restart": False,
                "write_restart_path": ic_base,
            },
            sources,
        )
        solve_steady_state(phys_steady)
        # Without a stored angular flux the initial angular flux is reconstructed from the flux
        # moments with the sources present at the switch to time-dependent mode, so the
        # steady-state source must still be set then; it is removed before stepping.
        phys_ic = make_problem({**common_options, "read_initial_condition_path": ic_base}, sources)
        ic_split = run_decay_steps(
            phys_ic,
            dt,
            dt,
            before_execute=lambda: phys_ic.SetVolumetricSources(clear_volumetric_sources=True),
        )

        # Transient restart: the writer's transient solver writes the restart after its step,
        # replacing the one its steady-state solve wrote.
        phys_writer = make_problem(
            {**common_options, "restart_writes_enabled": True, "write_restart_path": td_base},
            sources,
        )
        solve_steady_state(phys_writer)
        phys_writer.SetVolumetricSources(clear_volumetric_sources=True)
        phys_writer.SetTimeDependentMode()
        run_decay_steps(phys_writer, dt, dt)
        phys_reader = make_problem(
            {**common_options, "read_restart_path": td_base}, [], time_dependent=True
        )
        td_split = run_decay_steps(phys_reader, dt, 2.0 * dt)
    finally:
        for base in (ic_base, td_base):
            try:
                os.remove(f"{base}{rank}.restart.h5")
            except FileNotFoundError:
                pass

    return (
        relative_difference(ic_split, continuous_run(dt)),
        relative_difference(td_split, continuous_run(2.0 * dt)),
    )


def header_count_differences(grid, xs_maps, sources, quadrature, bases):
    # Each xs_map writes a steady-state restart without precursor data, and a problem on the
    # other xs_map reads it; read_restart_path checks the precursor count in the file header.
    options = {"verbose_inner_iterations": False, "verbose_outer_iterations": False}
    differences = []
    try:
        written = []
        for xs_map, base in zip(xs_maps, bases):
            phys = build_problem(
                grid,
                1,
                xs_map,
                sources,
                quadrature,
                {
                    **options,
                    "use_precursors": False,
                    "restart_writes_enabled": True,
                    "write_restart_path": base,
                },
            )
            solve_steady_state(phys)
            written.append(np.array(phys.GetPhiNewLocal(), copy=True))
        for writer, reader in ((0, 1), (1, 0)):
            if rank == 0:
                print(f"header_count: reading the restart of xs_map {writer} on xs_map {reader}")
            phys = build_problem(
                grid,
                1,
                xs_maps[reader],
                sources,
                quadrature,
                {**options, "use_precursors": True, "read_restart_path": bases[writer]},
            )
            SteadyStateSourceSolver(problem=phys).Initialize()
            phi_read = np.array(phys.GetPhiNewLocal(), copy=True)
            difference = global_max_abs([phi_read - written[writer]]) / global_max_abs(
                [written[writer]]
            )
            differences.append((writer, reader, difference))
    finally:
        for base in bases:
            try:
                os.remove(f"{base}{rank}.restart.h5")
            except FileNotFoundError:
                pass
    return differences


if __name__ == "__main__":
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")
    if "xs_case" not in globals() or xs_case not in XS_CASES:
        sys.exit(f"Set xs_case to one of {sorted(XS_CASES)}.")
    num_groups, assignments = XS_CASES[xs_case]

    nodes = [0.1 * i for i in range(11)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(infx=True, infy=True, zmin=0.0, zmax=0.5), 1, True
    )

    dt = 0.05
    source = VolumetricSource(
        block_ids=[0, 1], group_strength=[0.5] * num_groups, start_time=0.0, end_time=10.0
    )
    quadrature = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

    common_options = {
        "save_angular_flux": True,
        "use_precursors": True,
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": False,
    }

    # The restart files are written under out/, which git ignores, and each rank removes its own
    # files whether or not the run succeeds.
    restart_dir = os.path.join("out", f"transient_precursors_restart_two_blocks_{xs_case}")
    if rank == 0:
        shutil.rmtree(restart_dir, ignore_errors=True)
        os.makedirs(restart_dir)
    MPIBarrier()

    xs_maps = [load_xs_map(assignment) for assignment in assignments]
    max_diff = 0.0
    if xs_case == "header_count":
        bases = [os.path.join(restart_dir, f"steady_{index}_") for index in range(len(xs_maps))]
        for writer, reader, difference in header_count_differences(
            grid, xs_maps, [source], quadrature, bases
        ):
            if rank == 0:
                print(
                    f"header_count: xs_map {writer} restart read on xs_map {reader}: "
                    f"flux difference {difference:.6e}"
                )
            max_diff = max(max_diff, difference)
    else:
        for index, xs_map in enumerate(xs_maps):
            setup = (grid, num_groups, xs_map, [source], quadrature, common_options, dt)
            ic_diff, td_diff = restart_differences(
                setup,
                os.path.join(restart_dir, f"initial_condition_{index}_"),
                os.path.join(restart_dir, f"transient_{index}_"),
            )
            if rank == 0:
                print(
                    f"{xs_case} xs_map {index}: initial-condition restart difference "
                    f"{ic_diff:.6e}"
                )
                print(f"{xs_case} xs_map {index}: transient restart difference {td_diff:.6e}")
            max_diff = max(max_diff, ic_diff, td_diff)

    MPIBarrier()
    if rank == 0:
        shutil.rmtree(restart_dir, ignore_errors=True)
        print(f"PRECURSOR_RESTART_TWO_BLOCKS_DIFF={max_diff:.6e}")
