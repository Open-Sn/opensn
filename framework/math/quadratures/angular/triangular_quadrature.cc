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

  // For triangular quadrature, azimu_ang stores the maximum set of azimuthal angles
  // (at the equator), but each polar level may use a subset
  if (!azimuthal_per_polar.empty())
    azimu_ang = azimuthal_per_polar[Np / 2]; // Use equatorial level as reference

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
    new_omega.x = sin(qpoint.theta) * cos(qpoint.phi);
    new_omega.y = sin(qpoint.theta) * sin(qpoint.phi);
    new_omega.z = cos(qpoint.theta);

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

GLCTriangularQuadrature3DXYZ::GLCTriangularQuadrature3DXYZ(int Npolar,
                                                           int Nazimuthal,
                                                           int scattering_order,
                                                           bool verbose)
  : ProductQuadrature(3, scattering_order)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLCTriangularQuadrature3DXYZ: Npolar must be even.");

  if (Nazimuthal % 4 != 0)
    throw std::invalid_argument(
      "GLCTriangularQuadrature3DXYZ: Nazimuthal must be a multiple of 4.");

  // Check that we have enough azimuthal angles for the triangular reduction
  // The reduction is 4 azimuthal angles per level from equator (1 per octant)
  // Maximum reduction occurs at poles: 4 * (Npolar/2 - 1) levels from equator
  const int max_levels_from_equator = Npolar / 2 - 1;
  const int min_azimuthal_at_poles = Nazimuthal - 4 * max_levels_from_equator;
  if (min_azimuthal_at_poles < 4)
    throw std::invalid_argument(
      "GLCTriangularQuadrature3DXYZ: Nazimuthal too small for given Npolar. "
      "Need at least 4 azimuthal angles (1 per octant) at the poles. "
      "For Npolar=" +
      std::to_string(Npolar) + ", minimum Nazimuthal is " +
      std::to_string(4 * max_levels_from_equator + 4) + ".");

  GaussLegendreQuadrature gl_polar(Npolar);

  // Create polar angles
  std::vector<double> polar;
  for (int j = 0; j < Npolar; ++j)
    polar.emplace_back(M_PI - acos(gl_polar.qpoints[j][0]));

  // Calculate number of azimuthal angles per polar level
  // At equator: Nazimuthal angles
  // Each level away from equator: 4 less azimuthal angles (1 less per octant)
  // This maintains octant symmetry in the triangular pattern
  std::vector<int> num_azimuthal_per_level(Npolar);
  for (int j = 0; j < Npolar; ++j)
  {
    // Distance from equator in terms of polar levels
    // For even Npolar, the equatorial levels are at Npolar/2 - 1 and Npolar/2
    // Levels from equator: how many steps from the nearest equatorial level
    int levels_from_equator;
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

  for (int j = 0; j < Npolar; ++j)
  {
    const int Na_j = num_azimuthal_per_level[j];

    // Create azimuthal angles for this polar level (Chebyshev distribution)
    for (int i = 0; i < Na_j; ++i)
      azimuthal_per_polar[j].emplace_back(M_PI * (2 * (i + 1) - 1) / Na_j);

    // Create uniform weights for this polar level (all equal)
    // These will be normalized to sum to 1.0 in AssembleTriangularCosines
    for (int i = 0; i < Na_j; ++i)
      wts_per_polar[j].emplace_back(1.0);
  }

  // Initialize
  AssembleTriangularCosines(azimuthal_per_polar, polar, wts_per_polar, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();

  // Count total number of angles
  size_t total_angles = 0;
  for (int j = 0; j < Npolar; ++j)
    total_angles += num_azimuthal_per_level[j];

  log.Log() << "Using 3D XYZ Triangular quadrature with " << total_angles << " angles (" << Npolar
            << " polar levels, " << Nazimuthal << " max azimuthal)"
            << " and weight sum of " << std::fixed << std::setprecision(2) << weight_sum_;

  // Log the number of azimuthal angles at each polar level
  std::stringstream ss;
  ss << "Azimuthal angles per polar level: [";
  for (int j = 0; j < Npolar; ++j)
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

  // For triangular quadrature, azimu_ang stores the maximum set of azimuthal angles
  // (at the equator), but each polar level may use a subset
  if (!azimuthal_per_polar.empty())
    azimu_ang = azimuthal_per_polar[Np - 1]; // Use equatorial level (last in 2D) as reference

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
    new_omega.x = sin(qpoint.theta) * cos(qpoint.phi);
    new_omega.y = sin(qpoint.theta) * sin(qpoint.phi);
    new_omega.z = cos(qpoint.theta);

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

GLCTriangularQuadrature2DXY::GLCTriangularQuadrature2DXY(int Npolar,
                                                         int Nazimuthal,
                                                         int scattering_order,
                                                         bool verbose)
  : ProductQuadrature(2, scattering_order)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLCTriangularQuadrature2DXY: Npolar must be even.");

  if (Nazimuthal % 4 != 0)
    throw std::invalid_argument("GLCTriangularQuadrature2DXY: Nazimuthal must be a multiple of 4.");

  // For 2D, we only use the upper hemisphere (half of the polar levels)
  const int half = Npolar / 2;

  // Check that we have enough azimuthal angles for the triangular reduction
  // The reduction is 4 azimuthal angles per level from equator (1 per octant)
  // Maximum reduction occurs at the pole: (half - 1) levels from equator
  const int max_levels_from_equator = half - 1;
  const int min_azimuthal_at_pole = Nazimuthal - 4 * max_levels_from_equator;
  if (min_azimuthal_at_pole < 4)
    throw std::invalid_argument(
      "GLCTriangularQuadrature2DXY: Nazimuthal too small for given Npolar. "
      "Need at least 4 azimuthal angles (1 per octant) at the pole. "
      "For Npolar=" +
      std::to_string(Npolar) + ", minimum Nazimuthal is " +
      std::to_string(4 * max_levels_from_equator + 4) + ".");

  GaussLegendreQuadrature gl_polar(Npolar);

  // Create polar angles (only upper hemisphere - first half of GL points)
  std::vector<double> polar;
  for (int j = 0; j < half; ++j)
    polar.emplace_back(M_PI - acos(gl_polar.qpoints[j][0]));

  // Calculate number of azimuthal angles per polar level
  // At equator (j = half - 1): Nazimuthal angles
  // Each level toward the pole: 4 less azimuthal angles (1 less per octant)
  std::vector<int> num_azimuthal_per_level(half);
  for (int j = 0; j < half; ++j)
  {
    // Distance from equator: equator is at j = half - 1
    const int levels_from_equator = (half - 1) - j;
    // Reduce by 4 azimuthal angles per level (1 per octant) to maintain symmetry
    num_azimuthal_per_level[j] = Nazimuthal - 4 * levels_from_equator;
  }

  // Create azimuthal angles and weights per polar level
  std::vector<std::vector<double>> azimuthal_per_polar(half);
  std::vector<std::vector<double>> wts_per_polar(half);

  for (int j = 0; j < half; ++j)
  {
    const int Na_j = num_azimuthal_per_level[j];

    // Create azimuthal angles for this polar level (Chebyshev distribution)
    for (int i = 0; i < Na_j; ++i)
      azimuthal_per_polar[j].emplace_back(M_PI * (2 * (i + 1) - 1) / Na_j);

    // Create uniform weights for this polar level (all equal)
    // These will be normalized to sum to 1.0 in AssembleTriangularCosines
    for (int i = 0; i < Na_j; ++i)
      wts_per_polar[j].emplace_back(1.0);
  }

  // Initialize
  AssembleTriangularCosines(azimuthal_per_polar, polar, wts_per_polar, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();

  // Count total number of angles
  size_t total_angles = 0;
  for (int j = 0; j < half; ++j)
    total_angles += num_azimuthal_per_level[j];

  log.Log() << "Using 2D XY Triangular quadrature with " << total_angles << " angles (" << half
            << " polar levels, " << Nazimuthal << " max azimuthal)"
            << " and weight sum of " << std::fixed << std::setprecision(2) << weight_sum_;

  // Log the number of azimuthal angles at each polar level
  std::stringstream ss;
  ss << "Azimuthal angles per polar level: [";
  for (int j = 0; j < half; ++j)
  {
    if (j > 0)
      ss << ", ";
    ss << num_azimuthal_per_level[j];
  }
  ss << "]";
  log.Log() << ss.str();
}

} // namespace opensn
