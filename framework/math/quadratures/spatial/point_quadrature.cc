#include "framework/math/quadratures/spatial/point_quadrature.h"

namespace opensn
{

PointQuadrature::PointQuadrature()
  : SpatialQuadrature(QuadratureOrder::CONSTANT)
{
  qpoints_ = {{0.0, 0.0, 0.0}};
  weights_ = {1.0};
}

} // namespace opensn
