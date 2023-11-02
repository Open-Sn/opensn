#ifndef QUADRATURE_HEXAHEDRON_H
#define QUADRATURE_HEXAHEDRON_H

#include "quadrature.h"

namespace chi_math
{
  class QuadratureHexahedron;
}

//###################################################################
/**Quadrature set for tetrahedrons.*/
class chi_math::QuadratureHexahedron : public chi_math::Quadrature
{
public:
  //Constructor
  explicit
  QuadratureHexahedron(QuadratureOrder order);

};

#endif //QUADRATURE_HEXAHEDRON_H