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
public:
  ~ProductQuadrature() override = default;

  /// Return the abscissae index for the given polar and azimuthal angle indices.
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

  const std::vector<double>& GetPolarAngles() const { return polar_ang_; }

  const std::vector<double>& GetAzimuthalAngles() const { return azimu_ang_; }

  unsigned int GetNumPolarAngles() const override { return n_polar_; }

  unsigned int GetNumAzimuthalAngles() const override { return n_azimuthal_; }

protected:
  std::vector<double> polar_ang_;
  std::vector<double> azimu_ang_;

  ProductQuadrature(unsigned int dimension,
                    unsigned int scattering_order,
                    OperatorConstructionMethod method)
    : AngularQuadrature(
        AngularQuadratureType::PRODUCT_QUADRATURE, dimension, scattering_order, method),
      weight_sum_(0.0)
  {
  }

  /// Initializes the quadrature with custom angles and weights.
  void AssembleCosines(const std::vector<double>& azimuthal,
                       const std::vector<double>& polar,
                       const std::vector<double>& wts,
                       bool verbose);

  double weight_sum_;

  /// Linear indices of ordered directions mapped to polar level.
  std::map<unsigned int, std::vector<unsigned int>> map_directions_;

  /// Number of polar angles specified for this product quadrature.
  unsigned int n_polar_ = 0;

  /// Number of azimuthal angles specified for this product quadrature.
  unsigned int n_azimuthal_ = 0;
};

class GLProductQuadrature1DSlab : public ProductQuadrature
{
public:
  /// Constructor for 1D slab Gauss-Legendre product quadrature
  explicit GLProductQuadrature1DSlab(
    unsigned int Npolar,
    unsigned int scattering_order,
    bool verbose = false,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

  std::string GetName() const override { return "1D Slab Gauss-Legendre product"; }
};

class GLCProductQuadrature2DXY : public ProductQuadrature
{
public:
  /// Constructor for 2D XY Gauss-Legendre Chebyshev product quadrature
  explicit GLCProductQuadrature2DXY(
    unsigned int Npolar,
    unsigned int Nazimuthal,
    unsigned int scattering_order,
    bool verbose = false,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

  std::string GetName() const override { return "2D XY Gauss-Legendre/Chebyshev product"; }
};

class GLCProductQuadrature3DXYZ : public ProductQuadrature
{
public:
  /// Constructor for 3D XYZ Gauss-Legendre Chebyshev product quadrature
  explicit GLCProductQuadrature3DXYZ(
    unsigned int Npolar,
    unsigned int Nazimuthal,
    unsigned int scattering_order,
    bool verbose = false,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

  std::string GetName() const override { return "3D XYZ Gauss-Legendre/Chebyshev product"; }
};

} // namespace opensn
