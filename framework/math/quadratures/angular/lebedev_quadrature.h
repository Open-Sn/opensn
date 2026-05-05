// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include <vector>
#include <string>
#include <map>

namespace opensn
{

/// Lebedev quadrature for angular integration on the full unit sphere (3D XYZ).
///
/// Provides symmetric, high-order point sets on S^2. Point data is loaded from
/// predefined tabulated values keyed by order.
class LebedevQuadrature3DXYZ : public AngularQuadrature
{
public:
  /// Construct a 3D Lebedev quadrature.
  ///
  /// \param quadrature_order Order of the Lebedev quadrature set to load.
  /// \param scattering_order Scattering order for moment calculations.
  /// \param verbose Enable verbose output.
  /// \param method Harmonic operator construction method used to build the moment operators.
  LebedevQuadrature3DXYZ(unsigned int quadrature_order,
                         unsigned int scattering_order,
                         bool verbose = false,
                         OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

  std::string GetName() const override { return "3D XYZ Lebedev"; }

  unsigned int GetQuadratureOrder() const override { return quadrature_order_; }

private:
  /// Load quadrature points for the given order from predefined tabulated data.
  void LoadFromOrder(unsigned int quadrature_order, bool verbose = false);

  unsigned int quadrature_order_ = 0;
};

/// Lebedev quadrature restricted to the upper hemisphere (2D XY).
///
/// Only includes points with z >= 0. Points on the equator (z = 0) have their weights
/// halved since they are shared between the upper and lower hemispheres.
class LebedevQuadrature2DXY : public AngularQuadrature
{
public:
  /// Construct a 2D Lebedev quadrature (upper hemisphere only).
  ///
  /// \param quadrature_order Order of the Lebedev quadrature set to load.
  /// \param scattering_order Scattering order for moment calculations.
  /// \param verbose Enable verbose output.
  /// \param method Harmonic operator construction method used to build the moment operators.
  LebedevQuadrature2DXY(unsigned int quadrature_order,
                        unsigned int scattering_order,
                        bool verbose = false,
                        OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

  std::string GetName() const override { return "2D XY Lebedev"; }

  unsigned int GetQuadratureOrder() const override { return quadrature_order_; }

private:
  /// Load upper-hemisphere points for the given order, halving weights at the equator.
  void LoadFromOrder(unsigned int quadrature_order, bool verbose = false);

  unsigned int quadrature_order_ = 0;
};

} // namespace opensn
