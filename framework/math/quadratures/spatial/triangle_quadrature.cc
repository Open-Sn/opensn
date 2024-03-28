#include "framework/math/quadratures/spatial/triangle_quadrature.h"
#include <stdexcept>

namespace opensn
{

TriangleQuadrature::TriangleQuadrature(QuadratureOrder order) : SpatialQuadrature(order)
{
  switch (order)
  {
    case QuadratureOrder::CONSTANT:
    case QuadratureOrder::FIRST:
    {
      qpoints = {{3.333333333333332593e-01, 3.333333333333334259e-01, 0.0}};
      weights = {5.000000000000000000e-01};
      break;
    }
    case QuadratureOrder::SECOND:
    {
      qpoints = {{1.666666666666666019e-01, 1.666666666666667962e-01, 0.0},
                 {6.666666666666665186e-01, 1.666666666666666574e-01, 0.0},
                 {1.666666666666668517e-01, 6.666666666666667407e-01, 0.0}};
      weights = {1.666666666666666574e-01, 1.666666666666666574e-01, 1.666666666666666574e-01};
      break;
    }
    case QuadratureOrder::THIRD:
    case QuadratureOrder::FOURTH:
    {
      qpoints = {{4.459484909159648347e-01, 1.081030181680702890e-01, 0.0},
                 {4.459484909159648902e-01, 4.459484909159649457e-01, 0.0},
                 {1.081030181680701918e-01, 4.459484909159650567e-01, 0.0},
                 {9.157621350977065977e-02, 9.157621350977090957e-02, 0.0},
                 {8.168475729804584029e-01, 9.157621350977071528e-02, 0.0},
                 {9.157621350977096508e-02, 8.168475729804585139e-01, 0.0}};
      weights = {1.116907948390057498e-01,
                 1.116907948390057498e-01,
                 1.116907948390057498e-01,
                 5.497587182766094233e-02,
                 5.497587182766094233e-02,
                 5.497587182766094233e-02};
      break;
    }
    default:
    {
      throw std::invalid_argument(std::string(__FUNCTION__) +
                                  " Invalid triangular quadrature order");
    }
  }
}

} // namespace opensn
