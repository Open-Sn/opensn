#pragma once

#include "math/Quadratures/curvilinear_angular_quadrature.h"

#include "math/Quadratures/quadrature.h"

namespace chi_math
{
class SphericalAngularQuadrature;
}

/** Spherical product angular quadrature. */
class chi_math::SphericalAngularQuadrature : public chi_math::CurvilinearAngularQuadrature
{
  //  Methods
public:
  /** Effective constructor. Initialize with one-dimensional quadrature.
   *  If not already present in the quadrature, the method inserts
   *  the starting directions. */
  SphericalAngularQuadrature(const chi_math::Quadrature& quad_polar, const bool verbose = false);
  /** Default destructor. */
  virtual ~SphericalAngularQuadrature() = default;

  void MakeHarmonicIndices(unsigned int scattering_order, int dimension) override;

private:
  /** Initialize with one-dimensional quadrature. */
  void Initialize(const chi_math::Quadrature& quad_polar, const bool verbose = false);
  /** Initialize parametrizing factors of the spherical angular quadrature,
   *  starting from a fully initialized underlying product quadrature. */
  void InitializeParameters();
};
