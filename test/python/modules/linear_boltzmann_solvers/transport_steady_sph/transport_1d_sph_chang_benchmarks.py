#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Chang's six one-group spherical transport benchmarks in 1D spherical geometry.

Reference:
  B. Chang, "Six 1D Polar Transport Test Problems for the Deterministic and Monte-Carlo Method,"
  LLNL-TR-819680 (2021), https://www.osti.gov/servlets/purl/1769096

OpenSnExamples solves the same problems on 3D tetrahedral meshes. All six are purely absorbing
spheres of radius 1 with no scattering, so the exact solution is the uncollided flux:
  1. Interior Dirichlet:  sigma_t = 1, isotropic incident flux psi_0 = 1 per steradian.
  2. Radiating ball:      void, source q = 1 in r < 1/2, vacuum surface.
  3. Radiating shell:     void, source q = 1 in 1/3 < r < 2/3, vacuum surface.
  4. Pacman:              sigma_t = 1 in r < 1/2, void outside, incident psi_0 = 1.
  5. Barrier:             sigma_t = 1 in 1/3 < r < 2/3, void elsewhere, incident psi_0 = 1.
  6. Radiating absorber:  sigma_t = 1 and q = 1 in r < 1/2, void outside, vacuum surface.
Here q is the total isotropic emission density, so the source per steradian is q / (4 pi), and
the isotropic boundary group strength is 4 pi psi_0.

The reference scalar flux phi(r) = int psi(r, Omega) dOmega is evaluated independently of OpenSn
by tracing the exact uncollided angular flux along the backward ray through the piecewise-constant
shells and integrating over mu with Gauss-Legendre quadrature split at the interface tangency
cosines. Chang's expressions are mu-integrals of the same exponential attenuation, but most need
E2 or numerical quadrature (scipy); this single numpy-only evaluator covers all six problems.
It was cross-checked against Chang's expressions, as implemented with scipy in the OpenSnExamples
Problem_<k>.py inputs, at 97 radii in [0.013, 0.997] per problem: the maximum relative
difference is 7e-15 (problem 1), 5e-15 (2), 3e-15 (3), 2e-15 (4), 1e-15 (5), and 5e-15 (6).
The reference is then averaged over 12 equal-width shells with r^2 weighting (graded toward
the shell edges for the logarithmic singularities of problems 2 and 3); the shell averages
converge to better than 1e-8.

OpenSn solves each problem on 384 cells (32 per reference shell, with nodes at every material
interface) with an S256 Gauss-Legendre quadrature and computes the same shell averages with
FieldFunctionInterpolationVolume. The error is dominated by angular truncation: doubling the cells
changes it by less than 1%, and doubling the directions reduces it about fourfold. The check
SPH_CHANG_P<k>_ERR_S256 is the maximum relative shell-average error for problem k. It must match
the recorded value to within 1%: the result is bitwise identical on 1 to 4 ranks, and round-off
changes it by about 1e-8 of its value, so the margin admits platform variation while a change
in the discretization, quadrature, or source normalization fails.
"""

import math
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
    from pyopensn.aquad import GLProductQuadrature1DSpherical
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

# Problem: (shell radii, sigma_t per shell, q per shell, incident psi_0 per steradian)
PROBLEMS = {
    1: ([0.0, 1.0], [1.0], [0.0], 1.0),
    2: ([0.0, 0.5, 1.0], [0.0, 0.0], [1.0, 0.0], 0.0),
    3: ([0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0),
    4: ([0.0, 0.5, 1.0], [1.0, 0.0], [0.0, 0.0], 1.0),
    5: ([0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], 1.0),
    6: ([0.0, 0.5, 1.0], [1.0, 0.0], [1.0, 0.0], 0.0),
}
NUM_SHELLS = 12
CELLS_PER_SHELL = 32

_GL_X, _GL_W = np.polynomial.legendre.leggauss(16)


def exact_psi(r, mu, edges, sigma_t, q, psi_in):
    """Exact uncollided angular flux at radius r for the direction cosines mu (an array)."""
    p2 = r * r * (1.0 - mu * mu)
    s_out = r * mu + np.sqrt(np.maximum(edges[-1] ** 2 - p2, 0.0))
    # Distances along the backward ray to its crossings of every interior interface
    cuts = [np.zeros_like(mu), s_out]
    for rk in edges[1:-1]:
        d = np.sqrt(np.maximum(rk * rk - p2, 0.0))
        for s in (r * mu - d, r * mu + d):
            cuts.append(np.where((rk * rk > p2) & (s > 0.0) & (s < s_out), s, s_out))
    cuts = np.sort(np.array(cuts), axis=0)
    psi = np.zeros_like(mu)
    tau = np.zeros_like(mu)
    for s0, s1 in zip(cuts[:-1], cuts[1:]):
        length = s1 - s0
        s_mid = 0.5 * (s0 + s1)
        rho = np.sqrt(np.maximum(r * r - 2.0 * r * mu * s_mid + s_mid * s_mid, 0.0))
        shell = np.clip(np.searchsorted(edges, rho) - 1, 0, len(sigma_t) - 1)
        sig = np.asarray(sigma_t)[shell]
        optical = sig * length
        path = np.where(optical > 1.0e-12, -np.expm1(-optical) / np.where(sig > 0.0, sig, 1.0),
                        length)
        psi += np.asarray(q)[shell] / (4.0 * math.pi) * path * np.exp(-tau)
        tau += optical
    return psi + psi_in * np.exp(-tau)


def exact_phi(r, edges, sigma_t, q, psi_in):
    """Scalar flux 2 pi int psi dmu, split where the backward ray is tangent to an interface."""
    breaks = np.unique([-1.0, 0.0, 1.0] + [math.sqrt(1.0 - (rk / r) ** 2)
                                           for rk in edges[1:-1] if rk < r])
    t = 0.5 * (_GL_X + 1.0)
    total = 0.0
    for a, b in zip(breaks[:-1], breaks[1:]):
        # mu = a + (b - a) t^2 removes the sqrt(mu - a) behavior at a tangency cosine
        total += np.dot(_GL_W * (b - a) * t, exact_psi(r, a + (b - a) * t * t, edges, sigma_t,
                                                       q, psi_in))
    return 2.0 * math.pi * total


def exact_shell_average(r0, r1, edges, sigma_t, q, psi_in):
    """r^2-weighted average of the exact scalar flux over [r0, r1]."""
    grading = [0.0, 1.0 / 16.0, 1.0 / 8.0, 1.0 / 4.0, 0.5, 3.0 / 4.0, 7.0 / 8.0, 15.0 / 16.0, 1.0]
    total = 0.0
    for a, b in zip(grading[:-1], grading[1:]):
        for x, w in zip(_GL_X, _GL_W):
            r = r0 + (r1 - r0) * (a + 0.5 * (b - a) * (x + 1.0))
            total += 0.5 * (r1 - r0) * (b - a) * w * r * r * exact_phi(r, edges, sigma_t, q,
                                                                       psi_in)
    return 3.0 * total / (r1**3 - r0**3)


def zrange(z0, z1):
    return RPPLogicalVolume(infx=True, infy=True, zmin=float(z0), zmax=float(z1))


def solve_shell_averages(edges, sigma_t, q, psi_in, n_polar, shell_edges):
    nodes = np.linspace(0.0, 1.0, NUM_SHELLS * CELLS_PER_SHELL + 1).tolist()
    grid = OrthogonalMeshGenerator(node_sets=[nodes], coord_sys="spherical").Execute()
    xs_map = []
    sources = []
    for block_id, (r0, r1) in enumerate(zip(edges[:-1], edges[1:])):
        grid.SetBlockIDFromLogicalVolume(zrange(r0, r1), block_id, True)
        xs = MultiGroupXS()
        xs.CreateSimpleOneGroup(sigma_t[block_id], 0.0)
        xs_map.append({"block_ids": [block_id], "xs": xs})
        if q[block_id] > 0.0:
            sources.append(VolumetricSource(block_ids=[block_id], group_strength=[q[block_id]]))
    if psi_in > 0.0:
        bcs = [{"name": "zmax", "type": "isotropic", "group_strength": [4.0 * math.pi * psi_in]}]
    else:
        bcs = [{"name": "zmax", "type": "vacuum"}]

    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSpherical(n_polar=n_polar,
                                                                     scattering_order=0),
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 10,
            }
        ],
        xs_map=xs_map,
        volumetric_sources=sources,
        boundary_conditions=bcs,
        options={"verbose_inner_iterations": False},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    averages = []
    for r0, r1 in zip(shell_edges[:-1], shell_edges[1:]):
        avg = FieldFunctionInterpolationVolume()
        avg.SetOperationType("avg")
        avg.SetLogicalVolume(zrange(r0, r1))
        avg.AddFieldFunction(scalar_flux)
        avg.Execute()
        averages.append(avg.GetValue())
    return np.array(averages)


if __name__ == "__main__":
    if size != 2:
        sys.exit(f"Incorrect number of processors. Expected 2 but got {size}.")

    shell_edges = np.linspace(0.0, 1.0, NUM_SHELLS + 1)
    results = {}
    for k, (edges, sigma_t, q, psi_in) in PROBLEMS.items():
        numerical = solve_shell_averages(edges, sigma_t, q, psi_in, 256, shell_edges)
        if rank == 0:
            exact = np.array([exact_shell_average(r0, r1, np.array(edges), sigma_t, q, psi_in)
                              for r0, r1 in zip(shell_edges[:-1], shell_edges[1:])])
            results[k] = float(np.max(np.abs(numerical / exact - 1.0)))

    if rank == 0:
        for k, err in results.items():
            print(f"SPH_CHANG_P{k}_ERR_S256 {err:.6e}")
