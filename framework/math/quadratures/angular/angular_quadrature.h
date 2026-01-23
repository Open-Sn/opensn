// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include <memory>
#include <vector>

namespace opensn
{
struct QuadraturePointPhiTheta;

/// Angular quadrature type identifier.
enum class AngularQuadratureType
{
  ProductQuadrature = 1,
  SLDFEsq = 2,
  LebedevQuadrature = 3,
  TriangularQuadrature = 4,
};

/// Quadrature point in spherical coordinates.
struct QuadraturePointPhiTheta
{
  /// Azimuthal angle.
  double phi = 0.0;
  /// Polar angle.
  double theta = 0.0;

  QuadraturePointPhiTheta(const double phi, const double theta) : phi(phi), theta(theta) {}
};

/// Base class for angular quadratures used in discrete ordinates transport calculations.
class AngularQuadrature
{
public:
  /// Spherical harmonic indices for moment-to-direction mappings.
  struct HarmonicIndices
  {
    /// Degree of the spherical harmonic.
    unsigned int ell = 0;
    /// Order of the spherical harmonic.
    int m = 0;

    HarmonicIndices() = default;
    HarmonicIndices(unsigned int ell, int m) : ell(ell), m(m) {}

    bool operator==(const HarmonicIndices& other) const
    {
      return (ell == other.ell and m == other.m);
    }
  };

protected:
  explicit AngularQuadrature(AngularQuadratureType type,
                             unsigned int dimension,
                             unsigned int scattering_order)
    : type_(type), dimension_(dimension), scattering_order_(scattering_order)
  {
  }

  /// Discrete-to-moment operator matrix.
  std::vector<std::vector<double>> d2m_op_;
  /// Moment-to-discrete operator matrix.
  std::vector<std::vector<double>> m2d_op_;
  /// Mapping from moment index to spherical harmonic indices.
  std::vector<HarmonicIndices> m_to_ell_em_map_;
  /// Quadrature type identifier.
  AngularQuadratureType type_;
  /// Spatial dimension of the quadrature.
  unsigned int dimension_;
  /// Maximum scattering order for moment calculations.
  unsigned int scattering_order_;

  /// Populate the map of moment index to spherical harmonic indices.
  void MakeHarmonicIndices();

public:
  virtual ~AngularQuadrature() = default;

  /// Compute the discrete-to-moment operator.
  void BuildDiscreteToMomentOperator();

  /// Compute the moment-to-discrete operator.
  void BuildMomentToDiscreteOperator();

  /**
   * Return a reference to the precomputed discrete-to-moment operator.
   *
   * The operator is accessed as [d][m], where m is the moment index and d is the direction index.
   */
  std::vector<std::vector<double>> const& GetDiscreteToMomentOperator() const;

  /**
   * Return a reference to the precomputed moment-to-discrete operator.
   *
   * The operator is accessed as [d][m], where m is the moment index and d is the direction index.
   */
  std::vector<std::vector<double>> const& GetMomentToDiscreteOperator() const;

  /// Return a reference to the precomputed harmonic index map.
  const std::vector<HarmonicIndices>& GetMomentToHarmonicsIndexMap() const;

  unsigned int GetDimension() const { return dimension_; }

  unsigned int GetScatteringOrder() const { return scattering_order_; }

  unsigned int GetNumMoments() const { return m_to_ell_em_map_.size(); }

  AngularQuadratureType GetType() const { return type_; }

  /// Quadrature point abscissae in spherical coordinates.
  std::vector<QuadraturePointPhiTheta> abscissae;
  /// Quadrature weights.
  std::vector<double> weights;
  /// Quadrature point direction vectors in Cartesian coordinates.
  std::vector<Vector3> omegas;
};

} // namespace opensn
