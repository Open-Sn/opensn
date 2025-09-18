// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include <memory>
#include <vector>

namespace opensn
{
struct QuadraturePointPhiTheta;

enum class AngularQuadratureType
{
  ProductQuadrature = 1,
  SLDFEsq = 2,
  LebedevQuadrature = 3,
};

enum class OperatorConstructionMethod
{
  Standard = 0,          ///< Compute both operators directly
  GalerkinMethodOne = 1, ///< Compute M2D first, then D2M = inverse(M2D)
  GalerkinMethodTwo = 2  ///< Compute D2M first, then M2D = inverse(D2M)
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
  unsigned int dimension_;
  unsigned int scattering_order_;
  OperatorConstructionMethod construction_method_;

  explicit AngularQuadrature(
    AngularQuadratureType type,
    unsigned int dimension,
    unsigned int scattering_order,
    OperatorConstructionMethod method = OperatorConstructionMethod::Standard)
    : type_(type),
      dimension_(dimension),
      scattering_order_(scattering_order),
      construction_method_(method)
  {
  }

  /// Populates a map of moment m to the Spherical Harmonic indices required.
  void MakeHarmonicIndices();

  /// Helper method to invert a square matrix using PETSc
  std::vector<std::vector<double>> InvertMatrix(const std::vector<std::vector<double>>& matrix);

  /// Build D2M operator directly using spherical harmonics
  void BuildDiscreteToMomentOperatorStandard();

  /// Build M2D operator directly using spherical harmonics
  void BuildMomentToDiscreteOperatorStandard();

public:
  std::vector<QuadraturePointPhiTheta> abscissae;
  std::vector<double> weights;
  std::vector<Vector3> omegas;

  virtual ~AngularQuadrature() = default;

  /// Computes the discrete to moment operator based on construction method.
  void BuildDiscreteToMomentOperator();

  /// Computes the moment to discrete operator based on construction method.
  void BuildMomentToDiscreteOperator();

  /// Sets the operator construction method
  void SetOperatorConstructionMethod(OperatorConstructionMethod method)
  {
    construction_method_ = method;
  }

  /// Gets the current operator construction method
  OperatorConstructionMethod GetOperatorConstructionMethod() const { return construction_method_; }

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

  unsigned int GetDimension() const { return dimension_; }

  unsigned int GetScatteringOrder() const { return scattering_order_; }

  size_t GetNumMoments() const { return m_to_ell_em_map_.size(); }

  AngularQuadratureType GetType() const { return type_; }
};

} // namespace opensn
