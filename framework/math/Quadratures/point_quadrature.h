#pragma once

#include "framework/math/Quadratures/quadrature.h"

namespace chi_math
{

/**Quadrate for a single point. Helps generalize quadrature based integration
 * on 1D cell faces.*/
class PointQuadrature : public Quadrature
{
public:
  PointQuadrature();
};

} // namespace chi_math
