// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include <map>
#include <vector>

namespace opensn
{

/// Class for product quadratures
class ProductQuadrature : public AngularQuadrature
{
protected:
  double weight_sum_;

  /// Linear indices of ordered directions mapped to polar level.
  std::map<unsigned int, std::vector<unsigned int>> map_directions_;

  ProductQuadrature(int dimension, int scattering_order)
    : AngularQuadrature(AngularQuadratureType::ProductQuadrature, dimension, scattering_order),
      weight_sum_(0.0)
  {
  }

  /// Initializes the quadrature with custom angles and weights.
  void AssembleCosines(const std::vector<double>& azimuthal,
                       const std::vector<double>& polar,
                       const std::vector<double>& wts,
                       bool verbose);

public:
  std::vector<double> polar_ang;
  std::vector<double> azimu_ang;

  ~ProductQuadrature() override = default;

  /**
   * Obtains the abscissae index given the indices of the polar angle index and the azimuthal angle
   * index.
   */
  unsigned int GetAngleNum(const unsigned int polar_angle_index,
                           const unsigned int azimu_angle_index) const
  {
    return map_directions_.at(polar_angle_index)[azimu_angle_index];
  }

  /// Return constant reference to map_directions.
  const std::map<unsigned int, std::vector<unsigned int>>& GetDirectionMap() const
  {
    return map_directions_;
  }
};

class GLProductQuadrature1DSlab : public ProductQuadrature
{
public:
  /// Constructor for 1D slab Gauss-Legendre product quadrature
  explicit GLProductQuadrature1DSlab(int Npolar, int scattering_order, bool verbose = false);
};

class GLCProductQuadrature2DXY : public ProductQuadrature
{
public:
  /// Constructor for 2D XY Gauss-Legendre Chebyshev product quadrature
  explicit GLCProductQuadrature2DXY(int Npolar,
                                    int Nazimuthal,
                                    int scattering_order,
                                    bool verbose = false);
};

class GLCProductQuadrature3DXYZ : public ProductQuadrature
{
public:
  /// Constructor for 3D XYZ Gauss-Legendre Chebyshev product quadrature
  explicit GLCProductQuadrature3DXYZ(int Npolar,
                                     int Nazimuthal,
                                     int scattering_order,
                                     bool verbose = false);
};

} // namespace opensn
