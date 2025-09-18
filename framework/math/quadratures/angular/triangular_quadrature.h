// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include <map>
#include <vector>

namespace opensn
{

/// Base class for triangular quadratures.
///
/// Unlike product quadratures which have a fixed number of azimuthal angles per polar level,
/// triangular quadratures have a varying number that decreases as the polar angle moves away
/// from the equatorial plane.
class TriangularQuadrature : public AngularQuadrature
{
public:
  ~TriangularQuadrature() override = default;

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

  /// Polar angles for each polar level.
  std::vector<double> polar_ang;

protected:
  TriangularQuadrature(unsigned int dimension,
                       unsigned int scattering_order,
                       OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD)
    : AngularQuadrature(
        AngularQuadratureType::TRIANGULAR_QUADRATURE, dimension, scattering_order, method),
      weight_sum_(0.0)
  {
  }

  /// Sum of all quadrature weights.
  double weight_sum_;

  /// Azimuthal angles for each polar level.
  ///
  /// The number of azimuthal angles decreases as the polar angle moves away from the equator,
  /// in contrast to product quadratures which have a fixed count at every level.
  std::vector<std::vector<double>> azimuthal_per_polar_;

  /// Linear indices of ordered directions mapped to polar level.
  std::map<unsigned int, std::vector<unsigned int>> map_directions_;
};

/// Triangular GLC quadrature for 3D XYZ geometry.
///
/// Similar to GLCProductQuadrature3DXYZ, but with fewer azimuthal angles as the polar angle moves
/// away from the equatorial plane. The maximum number of azimuthal angles (at the equator) is
/// computed automatically as 2 * Npolar.
class GLCTriangularQuadrature3DXYZ : public TriangularQuadrature
{
public:
  /// Construct a 3D XYZ triangular Gauss-Legendre-Chebyshev quadrature.
  ///
  /// For each polar angle away from the equator there is 1 fewer azimuthal angle per octant.
  /// \param Npolar Number of polar angles.
  /// \param scattering_order Scattering order for moment calculations.
  /// \param verbose Enable verbose output.
  explicit GLCTriangularQuadrature3DXYZ(
    unsigned int Npolar,
    unsigned int scattering_order,
    bool verbose = false,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

private:
  /// Assemble quadrature points and weights for the varying azimuthal angles per polar level.
  void AssembleTriangularCosines(const std::vector<std::vector<double>>& azimuthal_per_polar,
                                 const std::vector<double>& polar,
                                 const std::vector<std::vector<double>>& wts_per_polar,
                                 bool verbose);
};

/// Triangular GLC quadrature for 2D XY geometry.
///
/// Similar to GLCProductQuadrature2DXY, but with fewer azimuthal angles as the polar angle moves
/// away from the equatorial plane. Only upper hemisphere points (z >= 0) are included. The maximum
/// number of azimuthal angles is computed automatically as 2 * Npolar.
class GLCTriangularQuadrature2DXY : public TriangularQuadrature
{
public:
  /// Construct a 2D XY triangular Gauss-Legendre-Chebyshev quadrature.
  ///
  /// For each polar angle away from the equator there is 1 fewer azimuthal angle per octant.
  /// \param Npolar Number of polar angles.
  /// \param scattering_order Scattering order for moment calculations.
  /// \param verbose Enable verbose output.
  explicit GLCTriangularQuadrature2DXY(
    unsigned int Npolar,
    unsigned int scattering_order,
    bool verbose = false,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

private:
  /// Assemble quadrature points and weights for the varying azimuthal angles per polar level.
  void AssembleTriangularCosines(const std::vector<std::vector<double>>& azimuthal_per_polar,
                                 const std::vector<double>& polar,
                                 const std::vector<std::vector<double>>& wts_per_polar,
                                 bool verbose);
};

} // namespace opensn
