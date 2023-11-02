#ifndef QUADRATURE_TETRAHEDRON_H
#define QUADRATURE_TETRAHEDRON_H

#include "quadrature.h"

namespace chi_math
{
class QuadratureTetrahedron;
}

//###################################################################
/**Quadrature set for tetrahedrons.*/
class chi_math::QuadratureTetrahedron : public chi_math::Quadrature
{
public:
  // Constructor
  explicit QuadratureTetrahedron(QuadratureOrder order);

  void KeastRule(const std::vector<std::vector<double>>& rule_data, const unsigned int n_pts);
};

#endif