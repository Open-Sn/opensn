#pragma once

#include "framework/math/quadratures/quadrature.h"

namespace chi_math
{
class QuadratureTriangle;
}

class chi_math::QuadratureTriangle : public chi_math::Quadrature
{
public:
  /**Initializes quadratures for use on triangles.*/
  explicit QuadratureTriangle(QuadratureOrder order);

  /*** The Dunavant rules are for triangles. This function takes
  permutation points and weights in a specific format as input and
  fills the _points and _weights vectors.*/
  void dunavant_rule(const double rule_data[][4], const unsigned int n_pts);

  void dunavant_rule2(const double* wts,
                      const double* a,
                      const double* b,
                      const unsigned int* permutation_ids,
                      const unsigned int n_wts);
};
