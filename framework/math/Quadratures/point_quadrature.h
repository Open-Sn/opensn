#ifndef CHITECH_POINT_QUADRATURE_H
#define CHITECH_POINT_QUADRATURE_H

#include "quadrature.h"

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

#endif // CHITECH_POINT_QUADRATURE_H
