#include "framework/math/quadratures/quadrature_quadrilateral.h"

#include "framework/math/quadratures/quadrature_gausslegendre.h"

namespace opensn
{

QuadratureQuadrilateral::QuadratureQuadrilateral(QuadratureOrder order) : Quadrature(order)
{
  QuadratureGaussLegendre legendre(order);

  legendre.SetRange({-1.0, 1.0});

  size_t N = legendre.qpoints_.size();

  qpoints_.resize(N * N);
  weights_.resize(N * N);

  unsigned int q = 0;

  for (unsigned int j = 0; j < N; j++)
    for (unsigned int i = 0; i < N; i++)
    {
      qpoints_[q](0) = legendre.qpoints_[i](0);
      qpoints_[q](1) = legendre.qpoints_[j](0);

      weights_[q] = legendre.weights_[i] * legendre.weights_[j];

      q++;
    }
}

} // namespace opensn
