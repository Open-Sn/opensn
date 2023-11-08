#pragma once

#include "framework/math/Quadratures/quadrature.h"

namespace chi_math
{
class QuadratureTetrahedron;
}

/**Quadrature set for tetrahedrons.*/
class chi_math::QuadratureTetrahedron : public chi_math::Quadrature
{
public:
  /**Initialzes a set of points for a quadrature integration over
   * the volume of a tetrahedron.*/
  explicit QuadratureTetrahedron(QuadratureOrder order);

  /** The Keast rules are for tets. This function takes permutation
  points and weights in a specific format as input and fills the
  _points and _weights vectors.*/
  void KeastRule(const std::vector<std::vector<double>>& rule_data, const unsigned int n_pts);
};
