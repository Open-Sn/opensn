#ifndef QUADRATURE_TRIANGLE_H
#define QUADRATURE_TRIANGLE_H

#include "quadrature.h"

namespace chi_math
{
class QuadratureTriangle;
}

class chi_math::QuadratureTriangle : public chi_math::Quadrature
{
public:
public:
  explicit QuadratureTriangle(QuadratureOrder order);

  void dunavant_rule(const double rule_data[][4], const unsigned int n_pts);

  void dunavant_rule2(const double* wts,
                      const double* a,
                      const double* b,
                      const unsigned int* permutation_ids,
                      const unsigned int n_wts);
};

#endif