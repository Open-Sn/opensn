// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/product_quadrature.h"
#include <map>
#include <vector>

namespace opensn
{

/// Triangular GLC quadrature for 3D XYZ geometry.
/// Similar to GLCProductQuadrature3DXYZ, but with fewer azimuthal angles
/// as the polar angle moves away from the equatorial plane.
class GLCTriangularQuadrature3DXYZ : public ProductQuadrature
{
public:
  /// Constructor for 3D XYZ Triangular Gauss-Legendre Chebyshev quadrature.
  /// For each polar angle away from the equator, there is 1 less azimuthal angle.
  explicit GLCTriangularQuadrature3DXYZ(int Npolar,
                                        int Nazimuthal,
                                        int scattering_order,
                                        bool verbose = false);

private:
  /// Assembles the quadrature points and weights for varying azimuthal angles per polar level.
  void AssembleTriangularCosines(const std::vector<std::vector<double>>& azimuthal_per_polar,
                                 const std::vector<double>& polar,
                                 const std::vector<std::vector<double>>& wts_per_polar,
                                 bool verbose);
};

/// Triangular GLC quadrature for 2D XY geometry.
/// Similar to GLCProductQuadrature2DXY, but with fewer azimuthal angles
/// as the polar angle moves away from the equatorial plane.
/// Only includes upper hemisphere points (z >= 0).
class GLCTriangularQuadrature2DXY : public ProductQuadrature
{
public:
  /// Constructor for 2D XY Triangular Gauss-Legendre Chebyshev quadrature.
  /// For each polar angle away from the equator, there is 1 less azimuthal angle per octant.
  explicit GLCTriangularQuadrature2DXY(int Npolar,
                                       int Nazimuthal,
                                       int scattering_order,
                                       bool verbose = false);

private:
  /// Assembles the quadrature points and weights for varying azimuthal angles per polar level.
  void AssembleTriangularCosines(const std::vector<std::vector<double>>& azimuthal_per_polar,
                                 const std::vector<double>& polar,
                                 const std::vector<std::vector<double>>& wts_per_polar,
                                 bool verbose);
};

} // namespace opensn
