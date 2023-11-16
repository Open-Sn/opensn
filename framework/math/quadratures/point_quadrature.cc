#include "framework/math/quadratures/point_quadrature.h"

namespace opensn
{

PointQuadrature::PointQuadrature() : Quadrature(QuadratureOrder::CONSTANT)
{
  qpoints_ = {{0.0, 0.0, 0.0}};
  weights_ = {1.0};
}

} // namespace opensn
