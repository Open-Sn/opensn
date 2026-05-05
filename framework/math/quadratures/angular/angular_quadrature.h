// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include "framework/data_types/ndarray.h"
#include <memory>
#include <string>
#include <vector>

namespace opensn
{
struct QuadraturePointPhiTheta;

/// Angular quadrature type identifier.
enum class AngularQuadratureType
{
  PRODUCT_QUADRATURE = 1,
  SLDFE_SQ = 2,
  LEBEDEV_QUADRATURE = 3,
  TRIANGULAR_QUADRATURE = 4,
};

/// Method used to construct the discrete-to-moment and moment-to-discrete operators.
enum class OperatorConstructionMethod
{
  STANDARD = 0,      ///< Standard weighted quadrature projection.
  GALERKIN_ONE = 1,  ///< Galerkin method 1: D2M = (M2D)^{-1}.
  GALERKIN_THREE = 3 ///< Galerkin method 3: orthogonalized approximate harmonics.
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
  explicit AngularQuadrature(
    AngularQuadratureType type,
    unsigned int dimension,
    unsigned int scattering_order,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD)
    : type_(type),
      dimension_(dimension),
      scattering_order_(scattering_order),
      construction_method_(method)
  {
  }

  /// Sets the operator construction method.
  void SetOperatorConstructionMethod(OperatorConstructionMethod method)
  {
    construction_method_ = method;
  }

  /// Populate the map of moment index to spherical harmonic indices.
  void MakeHarmonicIndices();
  /// Quadrature-specific order parameter used by harmonic-selection rules.
  virtual unsigned int GetQuadratureOrder() const { return 0; }
  /// Number of polar angles used by harmonic-selection rules.
  virtual unsigned int GetNumPolarAngles() const { return 0; }
  /// Number of azimuthal angles used by harmonic-selection rules.
  virtual unsigned int GetNumAzimuthalAngles() const { return 0; }

  /// Discrete-to-moment operator matrix.
  /// OpenSn storage is indexed [angle][moment].
  NDArray<double, 2> d2m_op_ = {};
  /// Moment-to-discrete operator matrix.
  /// OpenSn storage is indexed [angle][moment].
  NDArray<double, 2> m2d_op_ = {};
  /// Mapping from moment index to spherical harmonic indices.
  std::vector<HarmonicIndices> m_to_ell_em_map_ = {};
  /// Quadrature type identifier.
  AngularQuadratureType type_;
  /// Spatial dimension of the quadrature.
  unsigned int dimension_;
  /// Maximum scattering order for moment calculations.
  unsigned int scattering_order_;
  OperatorConstructionMethod construction_method_;
  /// Quadrature point abscissae in spherical coordinates.
  std::vector<QuadraturePointPhiTheta> abscissae_ = {};
  /// Quadrature weights.
  std::vector<double> weights_ = {};
  /// Quadrature point direction vectors in Cartesian coordinates.
  std::vector<Vector3> omegas_ = {};

public:
  virtual ~AngularQuadrature() = default;

  /// Compute the discrete-to-moment operator.
  void BuildDiscreteToMomentOperator();

  /// Compute the moment-to-discrete operator.
  void BuildMomentToDiscreteOperator();

  /// Gets the current operator construction method
  OperatorConstructionMethod GetOperatorConstructionMethod() const { return construction_method_; }

  /// Return a reference to the precomputed discrete-to-moment operator.
  /// The operator is accessed as (d,m), where d is the direction index and m is the moment index.
  NDArray<double, 2> const& GetDiscreteToMomentOperator() const;

  /// Return a reference to the precomputed moment-to-discrete operator.
  /// The operator is accessed as (d,m), where d is the direction index and m is the moment index.
  NDArray<double, 2> const& GetMomentToDiscreteOperator() const;

  /// Return a reference to the precomputed harmonic index map.
  const std::vector<HarmonicIndices>& GetMomentToHarmonicsIndexMap() const;

  unsigned int GetDimension() const { return dimension_; }

  unsigned int GetScatteringOrder() const { return scattering_order_; }

  unsigned int GetNumMoments() const { return m_to_ell_em_map_.size(); }

  size_t GetNumAngles() const { return omegas_.size(); }

  bool Empty() const { return GetNumAngles() == 0; }

  AngularQuadratureType GetType() const { return type_; }

  /// Return a user-facing quadrature name.
  virtual std::string GetName() const = 0;

  /// Return the sum of all quadrature weights.
  virtual double GetWeightSum() const;

  const std::vector<QuadraturePointPhiTheta>& GetAbscissae() const { return abscissae_; }

  const QuadraturePointPhiTheta& GetAbscissa(size_t angle_index) const
  {
    return abscissae_.at(angle_index);
  }

  const std::vector<double>& GetWeights() const { return weights_; }

  double GetWeight(size_t angle_index) const { return weights_.at(angle_index); }

  const std::vector<Vector3>& GetOmegas() const { return omegas_; }

  const Vector3& GetOmega(size_t angle_index) const { return omegas_.at(angle_index); }
};

} // namespace opensn
