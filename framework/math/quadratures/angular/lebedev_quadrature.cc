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
#include <cmath>
#include <stdexcept>

namespace opensn
{

LebedevQuadrature3DXYZ::LebedevQuadrature3DXYZ(int quadrature_order,
                                               unsigned int scattering_order,
                                               bool verbose)
  : AngularQuadrature(AngularQuadratureType::LebedevQuadrature, 3, scattering_order)
{
  LoadFromOrder(quadrature_order, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();
}

void
LebedevQuadrature3DXYZ::LoadFromOrder(int quadrature_order, bool verbose)
{
  abscissae.clear();
  weights.clear();
  omegas.clear();

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
    abscissae.push_back(qpoint);

    // Create the direction vector
    Vector3 omega{x / r, y / r, z / r};
    omegas.push_back(omega);

    // Store the weight
    weights.push_back(w);
    weight_sum += w;

    if (verbose)
    {
      char buf[200];
      snprintf(buf,
               200,
               "Varphi=%.2f Theta=%.2f Weight=%.3e\n",
               qpoint.phi * 180.0 / M_PI,
               qpoint.theta * 180.0 / M_PI,
               w);
      ostr << buf;
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
    for (auto& w : weights)
      w *= scale_factor;

    if (verbose)
      log.Log() << "Weights have been normalized to sum to 1.";
  }
}

} // namespace opensn
