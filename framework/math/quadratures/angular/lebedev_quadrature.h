// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include <vector>
#include <string>
#include <map>

namespace opensn
{

/**
 * @brief Implementation of Lebedev quadrature for angular integration on the surface of a sphere.
 *
 * Lebedev quadrature provides a set of points on the surface of a sphere that allows for
 * symmetric and efficient angular integration. This implementation reads point data from
 * files and constructs the quadrature set.
 */
class LebedevQuadrature3DXYZ : public AngularQuadrature
{
public:
  /**
   * @brief Constructor for Lebedev quadrature.
   *
   * @param quadrature_order The order of the Lebedev quadrature set to load
   * @param scattering_order Maximum scattering order
   * @param verbose Flag to enable verbose output
   * @param method Operator construction method
   */
  LebedevQuadrature3DXYZ(unsigned int quadrature_order,
                         unsigned int scattering_order,
                         bool verbose = false,
                         OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

private:
  /**
   * @brief Loads quadrature points for the specified order from predefined data.
   *
   * @param order The order to load
   * @param verbose Flag to enable verbose output
   */
  void LoadFromOrder(unsigned int quadrature_order, bool verbose = false);
};

} // namespace opensn
