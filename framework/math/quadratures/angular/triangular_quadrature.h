// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include <map>
#include <vector>

namespace opensn
{

/**
 * Base class for triangular quadratures.
 *
 * Unlike product quadratures which have a fixed number of azimuthal angles per polar level,
 * triangular quadratures have a varying number of azimuthal angles that decreases
 * as the polar angle moves away from the equatorial plane.
 */
class TriangularQuadrature : public AngularQuadrature
{
public:
  ~TriangularQuadrature() override = default;

  /// Obtain the abscissae index given the polar angle index and the azimuthal angle index.
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
        AngularQuadratureType::TriangularQuadrature, dimension, scattering_order, method),
      weight_sum_(0.0)
  {
  }

  /// Sum of all quadrature weights.
  double weight_sum_;

  /**
   * Azimuthal angles for each polar level.
   *
   * Unlike product quadratures which have a fixed number of azimuthal angles per polar level,
   * triangular quadratures have a varying number that decreases as the polar angle moves from
   * the equator.
   */
  std::vector<std::vector<double>> azimuthal_per_polar_;

  /// Linear indices of ordered directions mapped to polar level.
  std::map<unsigned int, std::vector<unsigned int>> map_directions_;
};

/**
 * Triangular GLC quadrature for 3D XYZ geometry.
 *
 * Similar to GLCProductQuadrature3DXYZ, but with fewer azimuthal angles as the polar angle moves
 * away from the equatorial plane. The number of azimuthal angles is maximized automatically based
 * on Npolar.
 */
class GLCTriangularQuadrature3DXYZ : public TriangularQuadrature
{
public:
  /**
   * Construct a 3D XYZ triangular Gauss-Legendre Chebyshev quadrature.
   *
   * For each polar angle away from the equator, there is 1 less azimuthal angle per octant.
   * The maximum number of azimuthal angles (at the equator) is computed as 2 * Npolar.
   *
   * \param Npolar Number of polar angles.
   * \param scattering_order Scattering order for moment calculations.
   * \param verbose Flag to enable verbose output.
   */
  explicit GLCTriangularQuadrature3DXYZ(
    unsigned int Npolar,
    unsigned int scattering_order,
    bool verbose = false,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

private:
  /**
   * Assemble the quadrature points and weights for varying azimuthal angles per polar level.
   *
   * \param azimuthal_per_polar Azimuthal angles for each polar level.
   * \param polar Polar angles.
   * \param wts_per_polar Weights for each polar level.
   * \param verbose Flag to enable verbose output.
   */
  void AssembleTriangularCosines(const std::vector<std::vector<double>>& azimuthal_per_polar,
                                 const std::vector<double>& polar,
                                 const std::vector<std::vector<double>>& wts_per_polar,
                                 bool verbose);
};

/**
 * Triangular GLC quadrature for 2D XY geometry.
 *
 * Similar to GLCProductQuadrature2DXY, but with fewer azimuthal angles as the polar angle moves
 * away from the equatorial plane. Only includes upper hemisphere points (z >= 0). The number of
 * azimuthal angles is maximized automatically based on Npolar.
 */
class GLCTriangularQuadrature2DXY : public TriangularQuadrature
{
public:
  /**
   * Construct a 2D XY triangular Gauss-Legendre Chebyshev quadrature.
   *
   * For each polar angle away from the equator, there is 1 less azimuthal angle per octant.
   * The maximum number of azimuthal angles (at the equator) is computed as 2 * Npolar.
   *
   * \param Npolar Number of polar angles.
   * \param scattering_order Scattering order for moment calculations.
   * \param verbose Flag to enable verbose output.
   */
  explicit GLCTriangularQuadrature2DXY(
    unsigned int Npolar,
    unsigned int scattering_order,
    bool verbose = false,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

private:
  /**
   * Assemble the quadrature points and weights for varying azimuthal angles per polar level.
   *
   * \param azimuthal_per_polar Azimuthal angles for each polar level.
   * \param polar Polar angles.
   * \param wts_per_polar Weights for each polar level.
   * \param verbose Flag to enable verbose output.
   */
  void AssembleTriangularCosines(const std::vector<std::vector<double>>& azimuthal_per_polar,
                                 const std::vector<double>& polar,
                                 const std::vector<std::vector<double>>& wts_per_polar,
                                 bool verbose);
};

} // namespace opensn
