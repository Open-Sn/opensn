#include "quadrature_wedge.h"

#include "quadrature_gausslegendre.h"
#include "quadrature_triangle.h"

namespace chi_math
{

QuadratureWedge::QuadratureWedge(QuadratureOrder order) : Quadrature(order)
{
  QuadratureGaussLegendre legendre(order);
  legendre.SetRange({-1.0, 1.0});
  QuadratureTriangle triangle(order);

  const size_t NL = legendre.qpoints_.size();
  const size_t NT = triangle.qpoints_.size();

  qpoints_.resize(NL * NT);
  weights_.resize(NL * NT);

  size_t q = 0;
  for (size_t i = 0; i < NL; ++i)
    for (size_t j = 0; j < NT; ++j)
    {
      qpoints_[q](0) = triangle.qpoints_[j][0];
      qpoints_[q](1) = triangle.qpoints_[j][1];
      qpoints_[q](2) = legendre.qpoints_[i][0];

      weights_[q] = legendre.weights_[i] * triangle.weights_[j];

      ++q;
    }
}

} // namespace chi_math