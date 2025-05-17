// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
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

LebedevQuadrature::LebedevQuadrature(int order, bool verbose)
  : AngularQuadrature(AngularQuadratureType::LebedevQuadrature, 3)
{
  LoadFromOrder(order, verbose);
}

void
LebedevQuadrature::LoadFromOrder(int order, bool verbose)
{
  abscissae.clear();
  weights.clear();
  omegas.clear();

  // Get points from LebedevOrders
  try
  {
    const auto& points = LebedevOrders::GetOrderPoints(order);

    std::stringstream ostr;
    double weight_sum = 0.0;
    int point_count = 0;

    for (const auto& point : points)
    {
      double x = point.x;
      double y = point.y;
      double z = point.z;
      double w = point.weight;

      // Calculate phi and theta from x, y, z
      double r = sqrt(x * x + y * y + z * z);
      double theta = acos(z / r);
      double phi = atan2(y, x);
      if (phi < 0.0)
        phi += 2.0 * M_PI;

      // Create the point
      QuadraturePointPhiTheta qpoint(phi, theta);
      abscissae.push_back(qpoint);

      // Create the direction vector
      Vector3 omega{x, y, z};
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

      point_count++;
    }

    // If no points found
    if (abscissae.empty())
      throw std::invalid_argument("LebedevQuadrature: No quadrature points found for order: " +
                                  std::to_string(order));

    if (verbose)
    {
      log.Log() << "Loaded " << point_count << " Lebedev quadrature points from order " << order;
      log.Log() << ostr.str() << "\n"
                << "Weight sum=" << weight_sum;
    }

    // Check weight sum (should be 1 for 3D)
    const double expected_sum = 1.0;
    if (fabs(weight_sum - expected_sum) > 1.0e-10)
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
  catch (const std::invalid_argument& e)
  {
    throw std::invalid_argument("LebedevQuadrature: Failed to load order " + std::to_string(order) +
                                ": " + e.what());
  }
}

} // namespace opensn