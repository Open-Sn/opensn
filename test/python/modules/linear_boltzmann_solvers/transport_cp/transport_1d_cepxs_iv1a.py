#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""CEPXS/ONELD Results Guide Problem IV.1.A.

Reference:
L. J. Lorence, Jr., J. E. Morel, and G. D. Valdez,
"Results Guide to CEPXS/ONELD: A One-Dimensional Coupled Electron-Photon
Discrete Ordinates Code Package," SAND89-2211, July 1990, p. 75.
https://www.osti.gov/servlets/purl/6758232

The problem is a normally incident 100 MeV electron beam on a 1 cm tungsten
slab. The 130-group library contains 50 electron, 30 photon, and 50 positron
groups with P15 scattering and full coupling. The CEPXS input used to generate
the checked-in BXSLIB file was::

    title
     W Brem
    cutoff 0.5
    energy 100.0
    legendre 15
    electrons
     linear 50
    photons
     log 30
    positrons
    electron-source
     full-coupling
    material w
    print-all

The BXSLIB SHA-256 is
436c30c079b3086e4bf541862f644a62b9be12a3b23ece6dec105cd67f0e0fa7.
The regression compares photon leakage per incident electron with the two
hemisphere-integrated discrete-ordinates values published in the report:
0.84 reverse and 8.20 forward. The tolerances account for the two-decimal
reference precision and for representing normal incidence with the most-forward
Gauss-Legendre ordinate.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.math import AngularFluxFunction
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.xs import MultiGroupXS


def make_normal_beam(quadrature):
    """Approximate a unit-current normal beam with the most-forward ordinate."""
    incoming = [
        (n, float(omega.z), float(quadrature.weights[n]))
        for n, omega in enumerate(quadrature.omegas)
        if float(omega.z) > 0.0 and float(quadrature.weights[n]) > 0.0
    ]
    if not incoming:
        raise RuntimeError("No incoming ordinates found for the zmin boundary")

    direction, mu, weight = max(incoming, key=lambda item: item[1])
    angular_flux = 1.0 / (mu * weight)

    def boundary_flux(group, angle):
        if group == 0 and angle == direction:
            return angular_flux
        return 0.0

    return AngularFluxFunction(boundary_flux)


def xs_path():
    return os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "../../../../assets/xs/W_50ge_30gp_50gpos_p15_CEPXS_standard.bxslib",
        )
    )


if __name__ == "__main__":
    num_cells = 50
    slab_thickness_cm = 1.0
    nodes = [i * slab_thickness_cm / num_cells for i in range(num_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromCEPXS(
        xs_path(),
        material_id=0,
        particle_order=["electron", "photon", "positron"],
    )
    if xs.num_groups != 130:
        raise RuntimeError(f"Expected 130 groups, found {xs.num_groups}")
    if xs.scattering_order != 15:
        raise RuntimeError(f"Expected P15 scattering, found P{xs.scattering_order}")

    quadrature = GLProductQuadrature1DSlab(
        n_polar=16,
        scattering_order=xs.scattering_order,
    )
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=xs.num_groups,
        groupsets=[
            {
                "groups_from_to": (0, xs.num_groups - 1),
                "angular_quadrature": quadrature,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 1000,
                "gmres_restart_interval": 300,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "zmin", "type": "arbitrary", "function": make_normal_beam(quadrature)},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={
            "field_function_prefix": "iv1a",
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    leakage = problem.ComputeLeakage(["zmin", "zmax"])
    photon_groups = slice(50, 80)
    reverse_photons = float(sum(leakage["zmin"][photon_groups]))
    forward_photons = float(sum(leakage["zmax"][photon_groups]))

    if rank == 0:
        print(f"IV1A_REVERSE_PHOTONS={reverse_photons:.12e}")
        print(f"IV1A_FORWARD_PHOTONS={forward_photons:.12e}")
