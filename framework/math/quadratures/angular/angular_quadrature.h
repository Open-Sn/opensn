// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/vector3.h"
#include <memory>
#include <vector>

namespace opensn
{
struct QuadraturePointPhiTheta;

enum class AngularQuadratureType
{
  ProductQuadrature = 1,
  SLDFESQ = 2
};

struct QuadraturePointPhiTheta
{
  double phi = 0.0;
  double theta = 0.0;

  QuadraturePointPhiTheta(const double phi, const double theta) : phi(phi), theta(theta) {}
};

class AngularQuadrature
{
public:
  struct HarmonicIndices
  {
    unsigned int ell = 0;
    int m = 0;

    HarmonicIndices() = default;
    HarmonicIndices(unsigned int ell, int m) : ell(ell), m(m) {}

    bool operator==(const HarmonicIndices& other) const
    {
      return (ell == other.ell and m == other.m);
    }
  };

protected:
  std::vector<std::vector<double>> d2m_op_;
  std::vector<std::vector<double>> m2d_op_;
  std::vector<HarmonicIndices> m_to_ell_em_map_;
  AngularQuadratureType type_;
  int dimension_;
  int scattering_order_;

  explicit AngularQuadrature(AngularQuadratureType type, int dimension, int scattering_order)
    : type_(type), dimension_(dimension), scattering_order_(scattering_order)
  {
  }

  /// Populates a map of moment m to the Spherical Harmonic indices required.
  void MakeHarmonicIndices();

public:
  std::vector<QuadraturePointPhiTheta> abscissae;
  std::vector<double> weights;
  std::vector<Vector3> omegas;

  virtual ~AngularQuadrature() = default;

  /// Computes the discrete to moment operator.
  void BuildDiscreteToMomentOperator();

  /// Computes the moment to discrete operator.
  void BuildMomentToDiscreteOperator();

  /**
   * Returns a reference to the precomputed d2m operator. The operator is accessed as [m][d], where
   * m is the moment index and d is the direction index.
   */
  std::vector<std::vector<double>> const& GetDiscreteToMomentOperator() const;

  /**
   * Returns a reference to the precomputed m2d operator. where m is the moment index and d is the
   * direction index.
   */
  std::vector<std::vector<double>> const& GetMomentToDiscreteOperator() const;

  /**
   * Returns a reference to the precomputed harmonic index map.
   */
  const std::vector<HarmonicIndices>& GetMomentToHarmonicsIndexMap() const;

  int GetDimension() { return dimension_; }

  int GetScatteringOrder() { return scattering_order_; }

  int GetNumMoments() { return m_to_ell_em_map_.size(); }

  AngularQuadratureType GetType() const { return type_; }
};

} // namespace opensn
