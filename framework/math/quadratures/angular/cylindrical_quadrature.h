#pragma once

#include "framework/math/quadratures/angular/curvilinear_quadrature.h"
#include "framework/math/quadratures/gauss_quadrature.h"

namespace opensn
{

class CylindricalQuadrature : public CurvilinearQuadrature
{
private:
  /** Initialize with one-dimensional quadratures: a polar quadrature and
   *  a possibly unique azimuthal quadrature for each polar level. */
  void Initialize(const GaussQuadrature& quad_polar,
                  const std::vector<GaussQuadrature>& quad_azimu_vec,
                  const bool verbose = false);

  /** Initialize parametrizing factors of the cylindrical angular quadrature,
   *  starting from a fully initialized underlying product quadrature. */
  void InitializeParameters();

public:
  /** Effective constructor. Initialize with one-dimensional quadratures:
   *  the azimuthal quadrature is applied at each polar level.
   *  If not already present in the azimuthal quadrature, the method inserts
   *  the azimuthal starting directions. */
  CylindricalQuadrature(const GaussQuadrature& quad_polar,
                        const GaussQuadrature& quad_azimu,
                        const bool verbose = false);

  /** Effective constructor. Initialize with one-dimensional quadratures:
   *  a possibly diverse azimuthal quadrature is applied at each polar level.
   *  If not already present in the azimuthal quadrature, the method inserts
   *  the azimuthal starting directions. */
  CylindricalQuadrature(const GaussQuadrature& quad_polar,
                        const std::vector<GaussQuadrature>& quad_azimu_vec,
                        const bool verbose = false);

  virtual ~CylindricalQuadrature() = default;

  void MakeHarmonicIndices(unsigned int scattering_order, int dimension) override;
};

} // namespace opensn
