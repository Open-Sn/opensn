#pragma once

#include "framework/math/quadratures/quadrature.h"

namespace opensn
{

/**Quadrate for a single point. Helps generalize quadrature based integration
 * on 1D cell faces.*/
class PointQuadrature : public Quadrature
{
public:
  PointQuadrature();
};

} // namespace opensn
