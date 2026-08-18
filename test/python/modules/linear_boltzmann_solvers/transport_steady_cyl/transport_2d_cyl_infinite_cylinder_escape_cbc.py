#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""RZ transport of an infinite cylinder converges to its first-flight solution.

An infinite cylinder of radius a = 1 cm holds a pure absorber, sigma_t = tau, with a uniform
isotropic source q = 1. With no scattering the scalar flux is the uncollided flux,

  phi(r) = (q / sigma_t) (1 - P(r)),
  P(r) = (1 / 4 pi) int exp(-sigma_t s(r, Omega)) dOmega,

where s(r, Omega) = d(r, varphi) / sin(theta) is the distance to the surface along Omega, and
d(r, varphi) = -r cos(varphi) + sqrt(a^2 - r^2 sin^2(varphi)) is that distance projected on the
cross section, with varphi the azimuth of Omega relative to the outward radial direction. The
flux falls from the axis to the surface (sigma_t phi = 0.726 at r = 0 and 0.496 at r = 0.9 for
tau = 1). The fraction of source particles that escape is P_esc = int_0^a 2 r P(r) dr / a^2.

Both references are evaluated here by Gauss-Legendre quadrature in mu = cos(theta), varphi and r
(96 points in each, which reproduces the 400-point values to 7e-7), independently of OpenSn.

Model: RZ mesh of 40 radial cells on [0, 1] and 2 axial cells on [0, 1], reflecting at the axis
and at both z planes (an infinite cylinder), vacuum at r = 1, GLC RZ quadratures. The escape
fraction is 1 - absorption / production from the balance table, and the flux is sampled at
r = 0, 0.25, 0.5, 0.75 and 0.9 on z = 0.5.

The error of these runs is set by the angular quadrature: with 40 radial cells the escape error
is 8e-4 for 8 polar and 16 azimuthal angles and 9e-5 (tau = 1) and 1.3e-4 (tau = 2) for 16 and
32. Refining both converges toward the reference (32 and 64 angles: 1.5e-5 with 40 cells and
4.5e-6 with 80). The tolerances are about 4 times the errors of the 16 x 32 runs.

Checks:
  RZ_ESCAPE_TAU<1|2>_REL_ERR: |P_esc / P_ref - 1| with 16 x 32 angles (tolerance 5e-4).
  RZ_FLUX_TAU<1|2>_MAX_REL_ERR: max over the sample points of |phi / phi_ref - 1| with 16 x 32
    angles (tolerance 4e-3; the largest error is at r = 0.9).
  RZ_ESCAPE_ANGULAR_CONVERGENCE: 1 when the tau = 1 escape error with 16 x 32 angles is at
    least 4 times smaller than with 8 x 16 angles (it is about 9 times smaller).
"""

import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.math import Vector3
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

NUM_RADIAL_CELLS = 40
SAMPLE_RADII = (0.0, 0.25, 0.5, 0.75, 0.9)
REFERENCE_POINTS = 96


def point_escape(tau, r, n=REFERENCE_POINTS):
    """P(r): probability that a particle born isotropically at radius r escapes uncollided."""
    mu, w_mu = np.polynomial.legendre.leggauss(n)
    x, w_x = np.polynomial.legendre.leggauss(2 * n)
    # varphi in (0, pi); the integrand is symmetric about varphi = 0.
    varphi = 0.5 * np.pi * (x + 1.0)
    w_varphi = 0.5 * np.pi * w_x
    sin_theta = np.sqrt(1.0 - mu**2)[:, None]
    d = (-r * np.cos(varphi) + np.sqrt(1.0 - (r * np.sin(varphi)) ** 2))[None, :]
    integral = np.sum(w_mu[:, None] * w_varphi[None, :] * np.exp(-tau * d / sin_theta))
    return float(integral / (2.0 * np.pi))


def escape_fraction(tau, n=REFERENCE_POINTS):
    """P_esc: the source-weighted average of P(r) over the cross section."""
    x, w = np.polynomial.legendre.leggauss(n)
    radii = 0.5 * (x + 1.0)
    return sum(0.5 * wi * 2.0 * ri * point_escape(tau, ri, n) for ri, wi in zip(radii, w))


def solve(tau, n_polar, n_azimuthal):
    """Returns the escape fraction and the scalar flux at the sample radii."""
    grid = OrthogonalMeshGenerator(
        node_sets=[[i / NUM_RADIAL_CELLS for i in range(NUM_RADIAL_CELLS + 1)], [0.0, 0.5, 1.0]],
        coord_sys="cylindrical",
    ).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(tau, 0.0)
    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DRZ(
                    n_polar=n_polar, n_azimuthal=n_azimuthal, scattering_order=0
                ),
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=[
            {"name": "rmin", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
            {"name": "rmax", "type": "vacuum"},
        ],
        options={"verbose_inner_iterations": False},
        sweep_type="CBC",
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    table = solver.ComputeBalanceTable()
    escape = 1.0 - table["absorption_rate"] / table["production_rate"]

    flux = problem.GetScalarFluxFieldFunction()[0]
    values = []
    for r in SAMPLE_RADII:
        interpolation = FieldFunctionInterpolationPoint()
        interpolation.SetPointOfInterest(Vector3(r, 0.5, 0.0))
        interpolation.AddFieldFunction(flux)
        interpolation.Execute()
        values.append(interpolation.GetPointValue())
    return escape, values


if __name__ == "__main__":
    results = {}
    escape_errors = {}
    for tau in (1.0, 2.0):
        escape_ref = escape_fraction(tau)
        flux_ref = [(1.0 - point_escape(tau, r)) / tau for r in SAMPLE_RADII]

        escape, flux = solve(tau, 16, 32)
        escape_errors[(tau, 16)] = abs(escape / escape_ref - 1.0)
        results[f"RZ_ESCAPE_TAU{int(tau)}_REL_ERR"] = escape_errors[(tau, 16)]
        results[f"RZ_FLUX_TAU{int(tau)}_MAX_REL_ERR"] = max(
            abs(value / ref - 1.0) for value, ref in zip(flux, flux_ref)
        )

    escape, _ = solve(1.0, 8, 16)
    coarse_error = abs(escape / escape_fraction(1.0) - 1.0)
    converged = int(coarse_error >= 4.0 * escape_errors[(1.0, 16)])

    if rank == 0:
        for key, value in results.items():
            print(f"{key}={value:.6e}")
        print(f"RZ_ESCAPE_ANGULAR_CONVERGENCE={converged}")
