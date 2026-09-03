#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, TransientSolver


def build_problem_options(grid, xs, source, quadrature, options, time_dependent=False):
    return {
        "mesh": grid,
        "num_groups": 1,
        "groupsets": [
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 200,
                "gmres_restart_interval": 30,
            }
        ],
        "xs_map": [{"block_ids": [0], "xs": xs}],
        "volumetric_sources": [source],
        "boundary_conditions": [
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        "options": options,
        "time_dependent": time_dependent,
        "sweep_type": "CBC",
    }


def evolve_transient(phys, dt, num_steps):
    solver = TransientSolver(
        problem=phys,
        dt=dt,
        theta=1.0,
        stop_time=dt * num_steps,
        initial_state="existing",
        verbose=False,
    )
    solver.Initialize()
    solver.Execute()

    return (
        np.array(phys.GetPhiNewLocal(), copy=True),
        [np.array(psi, copy=True) for psi in phys.GetPsi()],
    )


if __name__ == "__main__":
    meshgen = OrthogonalMeshGenerator(node_sets=[np.linspace(0.0, 1.0, 11).tolist()])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0, 1.0)
    source0 = VolumetricSource(block_ids=[0], group_strength=[2.0], start_time=0.0, end_time=10.0)
    quadrature = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

    dt = 0.01
    num_steps = 5

    common_options = {
        "save_angular_flux": True,
        "use_precursors": False,
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": False,
    }

    phys_continuous = DiscreteOrdinatesProblem(
        **build_problem_options(grid, xs, source0, quadrature, dict(common_options))
    )
    steady_continuous = SteadyStateSourceSolver(problem=phys_continuous)
    steady_continuous.Initialize()
    steady_continuous.Execute()
    phys_continuous.SetTimeDependentMode()
    phi_continuous, psi_continuous = evolve_transient(phys_continuous, dt, num_steps)

    for write_angular_flux in (False, True):
        for write_delayed_psi in (False, True):
            restart_base = f"steady_ic_cbc_af{int(write_angular_flux)}_dpsi{int(write_delayed_psi)}"
            steady_options = {
                **common_options,
                "restart_writes_enabled": True,
                "write_angular_flux_to_restart": write_angular_flux,
                "write_delayed_psi_to_restart": write_delayed_psi,
                "write_restart_path": restart_base,
            }
            phys_steady = DiscreteOrdinatesProblem(
                **build_problem_options(grid, xs, source0, quadrature, steady_options)
            )
            steady_solver = SteadyStateSourceSolver(problem=phys_steady)
            steady_solver.Initialize()
            steady_solver.Execute()

            transient_options = {
                **common_options,
                "read_initial_condition_path": restart_base,
            }
            phys_transient = DiscreteOrdinatesProblem(
                **build_problem_options(grid, xs, source0, quadrature, transient_options)
            )
            phi_split, psi_split = evolve_transient(phys_transient, dt, num_steps)

            local_ok = np.allclose(phi_continuous, phi_split, rtol=1.0e-10, atol=1.0e-12)
            phi_diff = np.max(np.abs(phi_continuous - phi_split))
            psi_diff = 0.0
            local_ok = local_ok and len(psi_continuous) == len(psi_split)
            for arr_continuous, arr_split in zip(psi_continuous, psi_split):
                psi_diff = max(psi_diff, np.max(np.abs(arr_continuous - arr_split)))
                local_ok = local_ok and np.allclose(
                    arr_continuous, arr_split, rtol=1.0e-10, atol=1.0e-12
                )

            global_ok = MPIAllReduce(float(local_ok))
            if global_ok != size:
                if rank == 0:
                    print(f"STEADY_RESTART_IC_WRITE_PSI {int(write_angular_flux)}")
                    print(f"STEADY_RESTART_IC_WRITE_DELAYED_PSI {int(write_delayed_psi)}")
                    print(f"STEADY_RESTART_IC_PHI_DIFF {phi_diff:.12e}")
                    print(f"STEADY_RESTART_IC_PSI_DIFF {psi_diff:.12e}")
                    sys.stdout.flush()
                raise ValueError("steady restart initial condition mismatch")

            restart_file = f"{restart_base}{rank}.restart.h5"
            if os.path.exists(restart_file):
                os.remove(restart_file)

    if rank == 0:
        print("STEADY_RESTART_IC_PASS 1")
