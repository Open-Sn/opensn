// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/lebedev_quadrature.h"
#include "framework/math/quadratures/angular/lebedev_orders.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace opensn
{

namespace
{

// Some Lebedev rules (for example orders 25 and 27) contain negative weights, which can produce
// negative angular fluxes and flux moments.
void
WarnIfNegativeWeights(const std::string& name,
                      unsigned int quadrature_order,
                      const std::vector<double>& weights)
{
  const auto num_negative =
    std::count_if(weights.begin(), weights.end(), [](double w) { return w < 0.0; });
  if (num_negative > 0)
    log.Log0Warning() << name << " of order " << quadrature_order << " has " << num_negative
                      << " negative weights. Angular fluxes and flux moments may become "
                         "negative; consider a different order.";
}

} // namespace

LebedevQuadrature3DXYZ::LebedevQuadrature3DXYZ(unsigned int quadrature_order,
                                               unsigned int scattering_order,
                                               bool verbose,
                                               OperatorConstructionMethod method)
  : AngularQuadrature(AngularQuadratureType::LEBEDEV_QUADRATURE, 3, scattering_order, method),
    quadrature_order_(quadrature_order)
{
  LoadFromOrder(quadrature_order, verbose);
  WarnIfNegativeWeights("LebedevQuadrature3DXYZ", quadrature_order, weights_);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();
}

void
LebedevQuadrature3DXYZ::LoadFromOrder(unsigned int quadrature_order, bool verbose)
{
  abscissae_.clear();
  weights_.clear();
  omegas_.clear();

  // Get points from LebedevOrders
  const auto& points = LebedevOrders::GetOrderPoints(quadrature_order);

  std::stringstream ostr;
  double weight_sum = 0.0;

  for (const auto& point : points)
  {
    const double x = point.x;
    const double y = point.y;
    const double z = point.z;
    const double w = point.weight;

    // Calculate phi and theta from x, y, z
    const double r = std::sqrt(x * x + y * y + z * z);
    const double theta = std::acos(z / r);
    double phi = std::atan2(y, x);
    if (phi < 0.0)
      phi += 2.0 * M_PI;

    // Create the point
    QuadraturePointPhiTheta qpoint(phi, theta);
    abscissae_.push_back(qpoint);

    // Create the direction vector
    Vector3 omega{x / r, y / r, z / r};
    omegas_.push_back(omega);

    // Store the weight
    weights_.push_back(w);
    weight_sum += w;

    if (verbose)
    {
      ostr << "Varphi=" << std::fixed << std::setprecision(2) << qpoint.phi * 180.0 / M_PI
           << " Theta=" << std::fixed << std::setprecision(2) << qpoint.theta * 180.0 / M_PI
           << " Weight=" << std::scientific << std::setprecision(3) << w << '\n';
    }
  }

  if (verbose)
  {
    log.Log() << "Loaded " << points.size() << " Lebedev quadrature points from quadrature order "
              << quadrature_order;
    log.Log() << ostr.str() << "\n"
              << "Weight sum=" << weight_sum;
  }

  // Check weight sum
  const double expected_sum = 1.0;
  if (std::fabs(weight_sum - expected_sum) > 1.0e-10)
  {
    if (verbose)
    {
      log.Log() << "Warning: Sum of weights differs from expected value 1.";
      log.Log() << "Expected: " << expected_sum << ", Actual: " << weight_sum;
    }

    // Normalize weights
    const double scale_factor = expected_sum / weight_sum;
    for (auto& w : weights_)
      w *= scale_factor;

    if (verbose)
      log.Log() << "Weights have been normalized to sum to 1.";
  }
}

LebedevQuadrature2DXY::LebedevQuadrature2DXY(unsigned int quadrature_order,
                                             unsigned int scattering_order,
                                             bool verbose,
                                             OperatorConstructionMethod method)
  : AngularQuadrature(AngularQuadratureType::LEBEDEV_QUADRATURE, 2, scattering_order, method),
    quadrature_order_(quadrature_order)
{
  LoadFromOrder(quadrature_order, verbose);
  WarnIfNegativeWeights("LebedevQuadrature2DXY", quadrature_order, weights_);

  // Every 2D Lebedev set includes the polar direction, which has no in-plane component.
  for (size_t n = 0; n < omegas_.size(); ++n)
    if (std::fabs(omegas_[n].x) < 1.0e-12 and std::fabs(omegas_[n].y) < 1.0e-12)
      log.Log0Warning() << "LebedevQuadrature2DXY: The set includes the polar direction (0, 0, 1) "
                        << "with weight " << weights_[n]
                        << ". It has no in-plane component, so in XY geometry it does not couple "
                           "to boundaries or neighboring cells and its angular flux is purely "
                           "local. In streaming-dominated regions this can bias the scalar flux by "
                           "up to the magnitude of that weight; the effect decreases with "
                           "quadrature order.";

  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();
}

void
LebedevQuadrature2DXY::LoadFromOrder(unsigned int quadrature_order, bool verbose)
{
  abscissae_.clear();
  weights_.clear();
  omegas_.clear();

  // Get points from LebedevOrders
  const auto& points = LebedevOrders::GetOrderPoints(quadrature_order);

  std::stringstream ostr;
  double weight_sum = 0.0;

  // Tolerance for determining if z is approximately zero (on the equator)
  const double z_tolerance = 1.0e-12;

  for (const auto& point : points)
  {
    const double x = point.x;
    const double y = point.y;
    const double z = point.z;
    double w = point.weight;

    // Skip points with z < 0 (lower hemisphere)
    if (z < -z_tolerance)
      continue;

    // For points on the equator (z ≈ 0), halve the weight
    if (std::fabs(z) <= z_tolerance)
      w *= 0.5;

    // Calculate phi and theta from x, y, z
    const double r = std::sqrt(x * x + y * y + z * z);
    const double theta = std::acos(z / r);
    double phi = std::atan2(y, x);
    if (phi < 0.0)
      phi += 2.0 * M_PI;

    // Create the point
    QuadraturePointPhiTheta qpoint(phi, theta);
    abscissae_.push_back(qpoint);

    // Create the direction vector
    Vector3 omega{x / r, y / r, z / r};
    omegas_.push_back(omega);

    // Store the weight (Doubled to renormalize from 0.5 to 1.0)
    weights_.push_back(2.0 * w);
    weight_sum += 2.0 * w;

    if (verbose)
    {
      ostr << "Varphi=" << std::fixed << std::setprecision(2) << qpoint.phi * 180.0 / M_PI
           << " Theta=" << std::fixed << std::setprecision(2) << qpoint.theta * 180.0 / M_PI
           << " Weight=" << std::scientific << std::setprecision(3) << w << '\n';
    }
  }

  if (verbose)
  {
    log.Log() << "Loaded " << GetNumAngles() << " Lebedev 2D quadrature points (upper hemisphere) "
              << "from quadrature order " << quadrature_order;
    log.Log() << ostr.str() << "\n"
              << "Weight sum=" << weight_sum;
  }

  // Check weight sum - should be 1.0 for upper hemisphere
  const double expected_sum = 1.0;
  if (std::fabs(weight_sum - expected_sum) > 1.0e-10)
  {
    if (verbose)
    {
      log.Log() << "Warning: Sum of weights differs from expected value 1.0.";
      log.Log() << "Expected: " << expected_sum << ", Actual: " << weight_sum;
    }
  }
}

} // namespace opensn
