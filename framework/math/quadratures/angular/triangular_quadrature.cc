// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/triangular_quadrature.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <cmath>
#include <sstream>
#include <iomanip>

namespace opensn
{

void
GLCTriangularQuadrature3DXYZ::AssembleTriangularCosines(
  const std::vector<std::vector<double>>& azimuthal_per_polar,
  const std::vector<double>& polar,
  const std::vector<std::vector<double>>& wts_per_polar,
  bool verbose)
{
  const size_t Np = polar.size();

  polar_ang = polar;
  azimuthal_per_polar_ = azimuthal_per_polar;

  if (verbose)
  {
    log.Log() << "Polar angles:";
    for (const auto& ang : polar_ang)
      log.Log() << ang;
  }

  // Create angle pairs with varying azimuthal angles per polar level
  map_directions_.clear();
  for (unsigned int j = 0; j < Np; ++j)
    map_directions_.emplace(j, std::vector<unsigned int>());

  abscissae.clear();
  weights.clear();
  weight_sum_ = 0.0;

  unsigned int direction_index = 0;
  for (unsigned int j = 0; j < Np; ++j)
  {
    const size_t Na_j = azimuthal_per_polar[j].size();

    if (verbose)
    {
      log.Log() << "Polar level " << j << " (theta=" << polar[j] * 180.0 / M_PI << " deg): " << Na_j
                << " azimuthal angles";
    }

    for (unsigned int i = 0; i < Na_j; ++i)
    {
      map_directions_[j].emplace_back(direction_index);

      const auto abscissa = QuadraturePointPhiTheta(azimuthal_per_polar[j][i], polar[j]);
      abscissae.emplace_back(abscissa);

      const double weight = wts_per_polar[j][i];
      weights.emplace_back(weight);
      weight_sum_ += weight;

      ++direction_index;
    }
  }

  // Create omega list
  omegas.clear();
  for (const auto& qpoint : abscissae)
  {
    Vector3 new_omega;
    new_omega.x = std::sin(qpoint.theta) * std::cos(qpoint.phi);
    new_omega.y = std::sin(qpoint.theta) * std::sin(qpoint.phi);
    new_omega.z = std::cos(qpoint.theta);

    omegas.emplace_back(new_omega);

    if (verbose)
      log.Log() << "Quadrature angle=" << new_omega.PrintStr();
  }

  // Normalize weights to 1.0
  for (auto& w : weights)
    w /= weight_sum_;

  // Compute and store sum of weights after normalization
  weight_sum_ = 0.0;
  for (auto& w : weights)
    weight_sum_ += w;
}

GLCTriangularQuadrature3DXYZ::GLCTriangularQuadrature3DXYZ(unsigned int Npolar,
                                                           unsigned int scattering_order,
                                                           bool verbose,
                                                           OperatorConstructionMethod method)
  : TriangularQuadrature(3, scattering_order, method)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLCTriangularQuadrature3DXYZ: Npolar must be even.");

  // Compute maximum number of azimuthal angles at the equator
  // For the triangular reduction pattern with 4 fewer azimuthal angles per level from equator,
  // this ensures exactly 4 azimuthal angles at the poles (1 per octant)
  const unsigned int Nazimuthal = 2 * Npolar;

  // Set parameters for Galerkin quadrature methods
  SetNumberOfPolar(Npolar);
  SetNumberOfAzimuthal(Nazimuthal);

  GaussLegendreQuadrature gl_polar(Npolar);

  // Create polar angles
  std::vector<double> polar;
  polar.reserve(Npolar);
  for (unsigned int j = 0; j < Npolar; ++j)
    polar.emplace_back(M_PI - std::acos(gl_polar.qpoints[j][0]));

  // Calculate number of azimuthal angles per polar level
  // At equator: Nazimuthal angles
  // Each level away from equator: 4 less azimuthal angles (1 less per octant)
  // This maintains octant symmetry in the triangular pattern
  std::vector<unsigned int> num_azimuthal_per_level(Npolar);
  for (unsigned int j = 0; j < Npolar; ++j)
  {
    // Distance from equator in terms of polar levels
    // For even Npolar, the equatorial levels are at Npolar/2 - 1 and Npolar/2
    // Levels from equator: how many steps from the nearest equatorial level
    unsigned int levels_from_equator = 0;
    if (j < Npolar / 2)
    {
      // Upper hemisphere: distance from equatorial level (Npolar/2 - 1)
      levels_from_equator = (Npolar / 2 - 1) - j;
    }
    else
    {
      // Lower hemisphere: distance from equatorial level (Npolar/2)
      levels_from_equator = j - (Npolar / 2);
    }
    // Reduce by 4 azimuthal angles per level (1 per octant) to maintain symmetry
    num_azimuthal_per_level[j] = Nazimuthal - 4 * levels_from_equator;
  }

  // Create azimuthal angles and weights per polar level
  std::vector<std::vector<double>> azimuthal_per_polar(Npolar);
  std::vector<std::vector<double>> wts_per_polar(Npolar);

  for (unsigned int j = 0; j < Npolar; ++j)
  {
    const unsigned int Na_j = num_azimuthal_per_level[j];

    // Create azimuthal angles for this polar level (Chebyshev distribution)
    for (unsigned int i = 0; i < Na_j; ++i)
      azimuthal_per_polar[j].emplace_back(M_PI * (2 * (i + 1) - 1) / Na_j);

    // Create weights: Gauss-Legendre weight (polar) * Gauss-Chebyshev weight (azimuthal)
    // Gauss-Chebyshev weight for integration over [0, 2π) is 2π/N
    const double gl_weight = gl_polar.weights[j];
    const double gc_weight = 2.0 * M_PI / Na_j;
    for (unsigned int i = 0; i < Na_j; ++i)
      wts_per_polar[j].emplace_back(gl_weight * gc_weight);
  }

  // Initialize
  AssembleTriangularCosines(azimuthal_per_polar, polar, wts_per_polar, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();

  // Count total number of angles
  size_t total_angles = 0;
  for (unsigned int j = 0; j < Npolar; ++j)
    total_angles += num_azimuthal_per_level[j];

  log.Log() << "Using 3D XYZ Triangular quadrature with " << total_angles << " angles (" << Npolar
            << " polar levels, " << Nazimuthal << " max azimuthal)"
            << " and weight sum of " << std::fixed << std::setprecision(2) << weight_sum_;

  // Log the number of azimuthal angles at each polar level
  std::stringstream ss;
  ss << "Azimuthal angles per polar level: [";
  for (unsigned int j = 0; j < Npolar; ++j)
  {
    if (j > 0)
      ss << ", ";
    ss << num_azimuthal_per_level[j];
  }
  ss << "]";
  log.Log() << ss.str();
}

void
GLCTriangularQuadrature2DXY::AssembleTriangularCosines(
  const std::vector<std::vector<double>>& azimuthal_per_polar,
  const std::vector<double>& polar,
  const std::vector<std::vector<double>>& wts_per_polar,
  bool verbose)
{
  const size_t Np = polar.size();

  polar_ang = polar;
  azimuthal_per_polar_ = azimuthal_per_polar;

  if (verbose)
  {
    log.Log() << "Polar angles:";
    for (const auto& ang : polar_ang)
      log.Log() << ang;
  }

  // Create angle pairs with varying azimuthal angles per polar level
  map_directions_.clear();
  for (unsigned int j = 0; j < Np; ++j)
    map_directions_.emplace(j, std::vector<unsigned int>());

  abscissae.clear();
  weights.clear();
  weight_sum_ = 0.0;

  unsigned int direction_index = 0;
  for (unsigned int j = 0; j < Np; ++j)
  {
    const size_t Na_j = azimuthal_per_polar[j].size();

    if (verbose)
    {
      log.Log() << "Polar level " << j << " (theta=" << polar[j] * 180.0 / M_PI << " deg): " << Na_j
                << " azimuthal angles";
    }

    for (unsigned int i = 0; i < Na_j; ++i)
    {
      map_directions_[j].emplace_back(direction_index);

      const auto abscissa = QuadraturePointPhiTheta(azimuthal_per_polar[j][i], polar[j]);
      abscissae.emplace_back(abscissa);

      const double weight = wts_per_polar[j][i];
      weights.emplace_back(weight);
      weight_sum_ += weight;

      ++direction_index;
    }
  }

  // Create omega list
  omegas.clear();
  for (const auto& qpoint : abscissae)
  {
    Vector3 new_omega;
    new_omega.x = std::sin(qpoint.theta) * std::cos(qpoint.phi);
    new_omega.y = std::sin(qpoint.theta) * std::sin(qpoint.phi);
    new_omega.z = std::cos(qpoint.theta);

    omegas.emplace_back(new_omega);

    if (verbose)
      log.Log() << "Quadrature angle=" << new_omega.PrintStr();
  }

  // Normalize weights to 1.0
  for (auto& w : weights)
    w /= weight_sum_;

  // Compute and store sum of weights after normalization
  weight_sum_ = 0.0;
  for (auto& w : weights)
    weight_sum_ += w;
}

GLCTriangularQuadrature2DXY::GLCTriangularQuadrature2DXY(unsigned int Npolar,
                                                         unsigned int scattering_order,
                                                         bool verbose,
                                                         OperatorConstructionMethod method)
  : TriangularQuadrature(2, scattering_order, method)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLCTriangularQuadrature2DXY: Npolar must be even.");

  // For 2D, we only use the upper hemisphere (half of the polar levels)
  const unsigned int half = Npolar / 2;

  // Compute maximum number of azimuthal angles at the equator
  // For the triangular reduction pattern with 4 fewer azimuthal angles per level from equator,
  // this ensures exactly 4 azimuthal angles at the pole (1 per quadrant)
  const unsigned int Nazimuthal = 2 * Npolar;

  // Set parameters for Galerkin quadrature methods
  SetNumberOfPolar(Npolar);
  SetNumberOfAzimuthal(Nazimuthal);

  GaussLegendreQuadrature gl_polar(Npolar);

  // Create polar angles (only upper hemisphere - first half of GL points)
  std::vector<double> polar;
  polar.reserve(half);
  for (unsigned int j = 0; j < half; ++j)
    polar.emplace_back(M_PI - acos(gl_polar.qpoints[j][0]));

  // Calculate number of azimuthal angles per polar level
  // At equator (j = half - 1): Nazimuthal angles
  // Each level toward the pole: 4 less azimuthal angles (1 less per quadrant)
  std::vector<unsigned int> num_azimuthal_per_level(half);
  for (unsigned int j = 0; j < half; ++j)
  {
    // Distance from equator: equator is at j = half - 1
    const unsigned int levels_from_equator = (half - 1) - j;
    // Reduce by 4 azimuthal angles per level (1 per quadrant) to maintain symmetry
    num_azimuthal_per_level[j] = Nazimuthal - 4 * levels_from_equator;
  }

  // Create azimuthal angles and weights per polar level
  std::vector<std::vector<double>> azimuthal_per_polar(half);
  std::vector<std::vector<double>> wts_per_polar(half);

  for (unsigned int j = 0; j < half; ++j)
  {
    const unsigned int Na_j = num_azimuthal_per_level[j];

    // Create azimuthal angles for this polar level (Chebyshev distribution)
    for (unsigned int i = 0; i < Na_j; ++i)
      azimuthal_per_polar[j].emplace_back(M_PI * (2 * (i + 1) - 1) / Na_j);

    // Create weights: Gauss-Legendre weight (polar) * Gauss-Chebyshev weight (azimuthal)
    // Gauss-Chebyshev weight for integration over [0, 2π) is 2π/N
    const double gl_weight = gl_polar.weights[j];
    const double gc_weight = 2.0 * M_PI / Na_j;
    for (unsigned int i = 0; i < Na_j; ++i)
      wts_per_polar[j].emplace_back(gl_weight * gc_weight);
  }

  // Initialize
  AssembleTriangularCosines(azimuthal_per_polar, polar, wts_per_polar, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();

  // Count total number of angles
  size_t total_angles = 0;
  for (unsigned int j = 0; j < half; ++j)
    total_angles += num_azimuthal_per_level[j];

  log.Log() << "Using 2D XY Triangular quadrature with " << total_angles << " angles (" << half
            << " polar levels, " << Nazimuthal << " max azimuthal)"
            << " and weight sum of " << std::fixed << std::setprecision(2) << weight_sum_;

  // Log the number of azimuthal angles at each polar level
  std::stringstream ss;
  ss << "Azimuthal angles per polar level: [";
  for (unsigned int j = 0; j < half; ++j)
  {
    if (j > 0)
      ss << ", ";
    ss << num_azimuthal_per_level[j];
  }
  ss << "]";
  log.Log() << ss.str();
}

} // namespace opensn
