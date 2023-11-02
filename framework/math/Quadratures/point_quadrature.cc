#include "point_quadrature.h"

namespace chi_math
{

PointQuadrature::PointQuadrature() : Quadrature(QuadratureOrder::CONSTANT)
{
  qpoints_ = {{0.0, 0.0, 0.0}};
  weights_ = {1.0};
}

} // namespace chi_math