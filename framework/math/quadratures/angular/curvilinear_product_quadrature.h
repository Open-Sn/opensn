// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"

namespace opensn
{

/**
 * Base class for curvilinear angular quadratures (product angular quadratures with additional
 * direction-dependent parameters).
 */
class CurvilinearProductQuadrature : public ProductQuadrature
{
protected:
  /// Factor to account for angular diamond differencing.
  std::vector<double> fac_diamond_difference_;

  /**
   * Factor to account for discretisation of the component of the streaming operator that contains
   * the angular derivative.
   */
  std::vector<double> fac_streaming_operator_;

  CurvilinearProductQuadrature(int dimension, int scattering_order)
    : ProductQuadrature(dimension, scattering_order)
  {
  }

public:
  const std::vector<double>& GetDiamondDifferenceFactor() const { return fac_diamond_difference_; }

  const std::vector<double>& GetStreamingOperatorFactor() const { return fac_streaming_operator_; }

  virtual ~CurvilinearProductQuadrature() = default;
};

class GLCProductQuadrature2DRZ : public CurvilinearProductQuadrature
{
private:
  /**
   * Initialize with one-dimensional quadratures: a polar quadrature and a possibly unique azimuthal
   * quadrature for each polar level.
   */
  void Initialize(const GaussQuadrature& quad_polar,
                  const std::vector<GaussQuadrature>& quad_azimu_vec,
                  const bool verbose = false);

  /**
   * Initialize parametrizing factors of the cylindrical angular quadrature, starting from a fully
   * initialized underlying product quadrature.
   */
  void InitializeParameters();

public:
  GLCProductQuadrature2DRZ(int Npolar, int Nazimuthal, int scattering_order, bool verbose = false);

  virtual ~GLCProductQuadrature2DRZ() = default;

  void MakeHarmonicIndices();
};

class GLProductQuadrature1DSpherical : public CurvilinearProductQuadrature
{
private:
  /// Initialize with one-dimensional quadrature.
  void Initialize(int Npolar, const bool verbose = false);

  /**
   * Initialize parametrizing factors of the spherical angular quadrature, starting from a fully
   * initialized underlying product quadrature.
   */
  void InitializeParameters();

public:
  GLProductQuadrature1DSpherical(int Npolar, int scattering_order, bool verbose = false);

  virtual ~GLProductQuadrature1DSpherical() = default;

  void MakeHarmonicIndices();
};

} // namespace opensn
