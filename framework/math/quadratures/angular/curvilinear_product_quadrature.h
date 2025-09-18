// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"

namespace opensn
{

/// Base class for curvilinear angular quadratures.
///
/// Extends product quadratures with direction-dependent parametrizing factors
/// needed for the curvilinear streaming operator.
class CurvilinearProductQuadrature : public ProductQuadrature
{
public:
  const std::vector<double>& GetDiamondDifferenceFactor() const { return fac_diamond_difference_; }

  const std::vector<double>& GetStreamingOperatorFactor() const { return fac_streaming_operator_; }

  ~CurvilinearProductQuadrature() override = default;

protected:
  CurvilinearProductQuadrature(unsigned int dimension,
                               unsigned int scattering_order,
                               OperatorConstructionMethod method)
    : ProductQuadrature(dimension, scattering_order, method)
  {
  }

  /// Factor to account for angular diamond differencing.
  std::vector<double> fac_diamond_difference_;

  /// Factor for discretization of the angular-derivative component of the streaming operator.
  std::vector<double> fac_streaming_operator_;
};

class GLCProductQuadrature2DRZ : public CurvilinearProductQuadrature
{
public:
  GLCProductQuadrature2DRZ(
    unsigned int Npolar,
    unsigned int Nazimuthal,
    unsigned int scattering_order,
    bool verbose = false,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

  ~GLCProductQuadrature2DRZ() override = default;

  void MakeHarmonicIndices();

private:
  /// Initialize from a polar quadrature and a per-polar-level azimuthal quadrature vector.
  void Initialize(const GaussQuadrature& quad_polar,
                  const std::vector<GaussQuadrature>& quad_azimu_vec,
                  bool verbose = false);

  /// Initialize cylindrical parametrizing factors from the underlying product quadrature.
  void InitializeParameters();
};

class GLProductQuadrature1DSpherical : public CurvilinearProductQuadrature
{
public:
  GLProductQuadrature1DSpherical(
    unsigned int Npolar,
    unsigned int scattering_order,
    bool verbose = false,
    OperatorConstructionMethod method = OperatorConstructionMethod::STANDARD);

  ~GLProductQuadrature1DSpherical() override = default;

  void MakeHarmonicIndices();

private:
  /// Initialize with a 1D polar quadrature.
  void Initialize(unsigned int Npolar, bool verbose = false);

  /// Initialize spherical parametrizing factors from the underlying product quadrature.
  void InitializeParameters();
};

} // namespace opensn
