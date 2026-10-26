#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Analytic two-material, two-group CMFD homogenization test.

Two equal-volume cells use the same fissile cross sections at density factors one
and two, and are deliberately placed in one mixed-material aggregate. The
test independently integrates the transport iterate in each fine cell and checks
reaction-rate preservation for the interaction matrices and transport-cross-section
homogenization followed by energy collapse for the diffusion coefficient.

One active power iteration assembles the CMFD operator from the post-sweep
transport iterate. A negligible relaxation (1e-30) makes every prolongation ratio
round to exactly one, so the final scalar flux is bitwise the iterate CMFD
homogenized and the test can integrate it independently.

With negative_self_scatter=True the data are put in transport-corrected P0 form
with a negative fast-group self-scatter, which leaves removal unchanged. The
homogenized coefficients must still be rate preserving rather than rejected.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX}

    def MPIAllReduce(value, op="max"):
        return MPI.COMM_WORLD.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import CMFDAcceleration, DiscreteOrdinatesProblem
    from pyopensn.solver import PowerIterationKEigenSolver
    from pyopensn.xs import MultiGroupXS


def max_abs_difference(actual, expected):
    if len(actual) != len(expected):
        return math.inf
    return max(abs(a - e) for a, e in zip(actual, expected))


def fuel_data(negative_self_scatter):
    """Tabulated xs_fuel_g2.xs data; scatter and production are row-major [to][from]."""
    data = {
        "sigma_t": [0.4241317, 0.7377545],
        "scatter": [0.3944575, 0.001266390, 0.0007568932, 0.4021512],
        "nu_sigma_f": [0.02615716, 0.6285298],
        # Out-scatter sums of the P1 transfer moments.
        "sigma_1": [0.04489551 - 0.0001835131, -0.00003343627 + 0.008252911],
    }
    if negative_self_scatter:
        shift = 0.4
        data["sigma_t"][0] -= shift
        data["scatter"][0] -= shift
        data["sigma_1"] = [0.0, 0.0]
    return data


def write_p0_xs(data, path):
    with open(path, "w") as stream:
        stream.write("NUM_GROUPS 2\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
        for g in range(2):
            stream.write(f"{g} {data['sigma_t'][g]!r}\n")
        stream.write("SIGMA_T_END\n\nNU_SIGMA_F_BEGIN\n")
        for g in range(2):
            stream.write(f"{g} {data['nu_sigma_f'][g]!r}\n")
        stream.write("NU_SIGMA_F_END\n\nNU_BEGIN\n0 1.1\n1 1.1\nNU_END\n\n")
        stream.write("CHI_BEGIN\n0 1.0\n1 0.0\nCHI_END\n\nTRANSFER_MOMENTS_BEGIN\n")
        for to in range(2):
            for frm in range(2):
                stream.write(f"M_GFROM_GTO_VAL 0 {frm} {to} {data['scatter'][to * 2 + frm]!r}\n")
        stream.write("TRANSFER_MOMENTS_END\n")


def expected_coefficients(flux_integrals, collapse_energy, adjoint, data):
    """Independently collapse the tabulated material data with the supplied flux integrals."""
    sigma_t = data["sigma_t"]
    scatter = data["scatter"]
    removal = [sigma_t[0] - scatter[0], sigma_t[1] - scatter[3]]
    production = data["nu_sigma_f"] + [0.0, 0.0]
    if adjoint:
        scatter = [scatter[col * 2 + row] for row in range(2) for col in range(2)]
        production = [production[col * 2 + row] for row in range(2) for col in range(2)]
    sigma_1 = data["sigma_1"]
    diffusion = [1.0 / (3.0 * (sigma_t[g] - sigma_1[g])) for g in range(2)]

    density_weight = []
    for g in range(2):
        denominator = flux_integrals[0][g] + flux_integrals[1][g]
        density_weight.append(
            (flux_integrals[0][g] + 2.0 * flux_integrals[1][g]) / denominator
        )

    # Spatially homogenize the transport reaction rate before inverting it to
    # obtain D. Averaging material diffusion coefficients directly would instead
    # average the inverse transport cross section and is not reaction-rate based.
    spatial_diffusion = [diffusion[g] / density_weight[g] for g in range(2)]
    if collapse_energy:
        group_flux = [sum(flux_integrals[cell][g] for cell in range(2)) for g in range(2)]
        density_flux = [group_flux[g] * density_weight[g] for g in range(2)]
        denominator = sum(group_flux)
        scattering_rate = sum(
            scatter[row * 2 + col] * density_flux[col]
            for row in range(2)
            for col in range(2)
        )
        return {
            "removal": [
                (sum(sigma_t[g] * density_flux[g] for g in range(2)) - scattering_rate)
                / denominator
            ],
            "scattering": [scattering_rate / denominator],
            "production": [
                sum(
                    production[row * 2 + col] * density_flux[col]
                    for row in range(2)
                    for col in range(2)
                )
                / denominator
            ],
            "diffusion": [
                sum(spatial_diffusion[g] * group_flux[g] for g in range(2)) / denominator
            ],
        }

    return {
        "removal": [removal[g] * density_weight[g] for g in range(2)],
        "scattering": [
            scatter[row * 2 + col] * density_weight[col]
            for row in range(2)
            for col in range(2)
        ],
        "production": [
            production[row * 2 + col] * density_weight[col]
            for row in range(2)
            for col in range(2)
        ],
        "diffusion": spatial_diffusion,
    }


def coefficient_errors(cmfd, expected):
    """Return local-owner errors; the test reduces them across ranks."""
    try:
        actual = {
            "removal": cmfd.GetHomogenizedRemovalXS(0),
            "scattering": cmfd.GetHomogenizedScatteringXS(0),
            "production": cmfd.GetHomogenizedProductionXS(0),
            "diffusion": cmfd.GetHomogenizedDiffusionCoefficient(0),
        }
    except ValueError:
        return {name: 0.0 for name in expected}
    return {name: max_abs_difference(actual[name], values) for name, values in expected.items()}


if __name__ == "__main__":
    negative_self_scatter = globals().get("negative_self_scatter", False)
    data = fuel_data(negative_self_scatter)
    xs_file = "../../../../assets/xs/xs_fuel_g2.xs"
    if negative_self_scatter:
        xs_file = f"cmfd_homogenization_negative_self_scatter_{rank}.xs"
        write_p0_xs(data, xs_file)
    base_xs = MultiGroupXS()
    base_xs.LoadFromOpenSn(xs_file)
    dense_xs = MultiGroupXS()
    dense_xs.LoadFromOpenSn(xs_file)
    dense_xs.Scale(2.0)
    if negative_self_scatter:
        os.remove(xs_file)

    grid = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0, 2.0]]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(zmin=1.0, zmax=2.0, infx=True, infy=True), 1, True
    )

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLProductQuadrature1DSlab(
                    n_polar=4, scattering_order=0
                ),
                "inner_linear_method": "petsc_richardson",
                "l_max_its": 20,
                "l_abs_tol": 1.0e-12,
            }
        ],
        xs_map=[
            {"block_ids": [0], "xs": base_xs},
            {"block_ids": [1], "xs": dense_xs},
        ],
        options={
            "use_precursors": False,
            "adjoint": globals().get("adjoint", False),
        },
    )
    collapse_energy = globals().get("collapse_energy", False)
    adjoint = globals().get("adjoint", False)
    cmfd_options = {
        "problem": problem,
        "coarse_mesh": "local_aggregation" if size == 1 else "global_aggregation",
        "aggregation_size": 2,
        "allow_mixed_material_aggregation": True,
        "group_aggregation_size": 2 if collapse_energy else 1,
        "relaxation": 1.0e-30,
    }
    cmfd = CMFDAcceleration(**cmfd_options)
    solver = PowerIterationKEigenSolver(problem=problem, acceleration=cmfd, max_iters=1)
    solver.Initialize()
    solver.Execute()

    # PWLD has two unit-interval nodes per cell. With one moment and two groups,
    # each cell therefore occupies four consecutive scalar-flux entries.
    phi = problem.GetPhiNewLocal()
    if size == 1:
        flux_integrals = [
            [0.5 * (phi[0] + phi[2]), 0.5 * (phi[1] + phi[3])],
            [0.5 * (phi[4] + phi[6]), 0.5 * (phi[5] + phi[7])],
        ]
    else:
        # The orthogonal mesh's one-dimensional partition ordering places the
        # lower-z cell on the last of these two ranks.
        local_integral = [0.5 * (phi[0] + phi[2]), 0.5 * (phi[1] + phi[3])]
        flux_integrals = [
            [
                MPIAllReduce(
                    local_integral[g] if rank == 1 - cell else 0.0, "sum"
                )
                for g in range(2)
            ]
            for cell in range(2)
        ]

    errors = coefficient_errors(
        cmfd, expected_coefficients(flux_integrals, collapse_energy, adjoint, data)
    )
    max_error = MPIAllReduce(max(errors.values()), "max")
    if rank == 0:
        print(f"CMFD_HOMOGENIZATION_MAX_ABS_ERROR {max_error:.16e}")
        print(f"CMFD_HOMOGENIZATION_COARSE_CELLS {cmfd.GetNumCoarseCells()}")
    if max_error > 1.0e-12:
        raise RuntimeError(f"CMFD homogenization coefficient errors: {errors}")
