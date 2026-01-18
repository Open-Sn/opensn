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
   * @param order The order of the Lebedev quadrature set to load
   * @param verbose Flag to enable verbose output
   */
  LebedevQuadrature3DXYZ(unsigned int quadrature_order,
                         unsigned int scattering_order,
                         bool verbose = false);

private:
  /**
   * @brief Loads quadrature points for the specified order from predefined data.
   *
   * @param order The order to load
   * @param verbose Flag to enable verbose output
   */
  void LoadFromOrder(unsigned int quadrature_order, bool verbose = false);
};

/**
 * @brief Implementation of 2D Lebedev quadrature for angular integration on the upper hemisphere.
 *
 * This is a 2D version of the Lebedev quadrature that only includes points with z >= 0
 * (i.e., polar angles θ <= π/2). Points on the equator (z = 0) have their weights halved
 * since they are shared between the upper and lower hemispheres.
 */
class LebedevQuadrature2DXY : public AngularQuadrature
{
public:
  /**
   * @brief Constructor for 2D Lebedev quadrature.
   *
   * @param quadrature_order The order of the Lebedev quadrature set to load
   * @param scattering_order The scattering order for moment calculations
   * @param verbose Flag to enable verbose output
   */
  LebedevQuadrature2DXY(unsigned int quadrature_order,
                        unsigned int scattering_order,
                        bool verbose = false);

private:
  /**
   * @brief Loads quadrature points for the specified order from predefined data,
   *        keeping only upper hemisphere points (z >= 0) and halving weights for z = 0.
   *
   * @param quadrature_order The order to load
   * @param verbose Flag to enable verbose output
   */
  void LoadFromOrder(unsigned int quadrature_order, bool verbose = false);
};

} // namespace opensn
