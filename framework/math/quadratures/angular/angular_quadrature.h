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

class AngularQuadrature : public std::enable_shared_from_this<AngularQuadrature>
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
  bool d2m_op_built_ = false;
  bool m2d_op_built_ = false;
  AngularQuadratureType type_;
  int dimension_;

  explicit AngularQuadrature(AngularQuadratureType type, int dimension)
    : type_(type), dimension_(dimension)
  {
  }

  /// Populates a map of moment m to the Spherical Harmonic indices required.
  virtual void MakeHarmonicIndices(unsigned int scattering_order);

public:
  std::vector<QuadraturePointPhiTheta> abscissae;
  std::vector<double> weights;
  std::vector<Vector3> omegas;

  virtual ~AngularQuadrature() = default;

  /// Computes the discrete to moment operator.
  virtual void BuildDiscreteToMomentOperator(unsigned int scattering_order);

  /// Computes the moment to discrete operator.
  virtual void BuildMomentToDiscreteOperator(unsigned int scattering_order);

  /**
   * Returns a reference to the precomputed d2m operator. This will throw a std::logic_error if the
   * operator has not been built yet. The operator is accessed as [m][d], where m is the moment
   * index and d is the direction index.
   */
  std::vector<std::vector<double>> const& GetDiscreteToMomentOperator() const;

  /**
   * Returns a reference to the precomputed m2d operator. This will throw a std::logic_error if the
   * operator has not been built yet. The operator is accessed as [m][d], where m is the moment
   * index and d is the direction index.
   */
  std::vector<std::vector<double>> const& GetMomentToDiscreteOperator() const;

  /**
   * Returns a reference to the precomputed harmonic index map. This will throw a std::logic_error
   * if the map has not been built yet.
   */
  const std::vector<HarmonicIndices>& GetMomentToHarmonicsIndexMap() const;

  int GetDimension() { return dimension_; }

  AngularQuadratureType GetType() const { return type_; }
};

} // namespace opensn
