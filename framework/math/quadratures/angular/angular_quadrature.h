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
  unsigned int quadrature_order_ = 0;
  unsigned int n_polar_ = 0;
  unsigned int n_azimuthal_ = 0;

  /// Populate the map of moment index to spherical harmonic indices.
  void MakeHarmonicIndices();

public:
  virtual ~AngularQuadrature() = default;

  /// Compute the discrete-to-moment operator.
  void BuildDiscreteToMomentOperator();

  /// Compute the moment-to-discrete operator.
  void BuildMomentToDiscreteOperator();

  /// Sets the operator construction method
  void SetOperatorConstructionMethod(OperatorConstructionMethod method)
  {
    construction_method_ = method;
  }

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

  AngularQuadratureType GetType() const { return type_; }

  /// Return a user-facing quadrature name.
  virtual std::string GetName() const = 0;

  /// Return the sum of all quadrature weights.
  virtual double GetWeightSum() const;

  /// Set the quadrature-specific order parameter.
  /// For Lebedev: Lebedev order {3, 5, 7, ...}; for SLDFESQ: uniform refinement level.
  void SetQuadratureOrder(unsigned int order) { quadrature_order_ = order; }

  /// Set the number of positive-mu polar points (product quadrature types).
  void SetNumberOfPolar(unsigned int num_polar) { n_polar_ = num_polar; }

  /// Set the total number of azimuthal points (product quadrature types).
  void SetNumberOfAzimuthal(unsigned int num_azimu) { n_azimuthal_ = num_azimu; }

  /// Quadrature point abscissae in spherical coordinates.
  std::vector<QuadraturePointPhiTheta> abscissae = {};
  /// Quadrature weights.
  std::vector<double> weights = {};
  /// Quadrature point direction vectors in Cartesian coordinates.
  std::vector<Vector3> omegas = {};
};

} // namespace opensn
