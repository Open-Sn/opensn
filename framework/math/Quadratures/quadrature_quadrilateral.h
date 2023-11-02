#ifndef QUADRATURE_QUADRILATERAL_H
#define QUADRATURE_QUADRILATERAL_H

#include "quadrature.h"

namespace chi_math
{
class QuadratureQuadrilateral;
}

//###################################################################
/**Quadrature set for quadrilaterals.*/
class chi_math::QuadratureQuadrilateral : public chi_math::Quadrature
{
public:
  // Constructor
  explicit QuadratureQuadrilateral(QuadratureOrder order);
};

#endif // QUADRATURE_QUADRILATERAL_H