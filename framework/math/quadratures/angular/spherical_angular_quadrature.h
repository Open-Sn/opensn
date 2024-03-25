#pragma once

#include "framework/math/quadratures/angular/curvilinear_angular_quadrature.h"
#include "framework/math/quadratures/spatial/spatial_quadrature.h"

namespace opensn
{

/** Spherical product angular quadrature. */
class SphericalAngularQuadrature : public CurvilinearAngularQuadrature
{
  //  Methods
public:
  /** Effective constructor. Initialize with one-dimensional quadrature.
   *  If not already present in the quadrature, the method inserts
   *  the starting directions. */
  SphericalAngularQuadrature(const Quadrature& quad_polar, const bool verbose = false);
  /** Default destructor. */
  virtual ~SphericalAngularQuadrature() = default;

  void MakeHarmonicIndices(unsigned int scattering_order, int dimension) override;

private:
  /** Initialize with one-dimensional quadrature. */
  void Initialize(const Quadrature& quad_polar, const bool verbose = false);
  /** Initialize parametrizing factors of the spherical angular quadrature,
   *  starting from a fully initialized underlying product quadrature. */
  void InitializeParameters();
};

} // namespace opensn
