#include "framework/math/Quadratures/Conical/quadrature_conical.h"

#include "framework/math/Quadratures/quadrature_gausslegendre.h"
#include "framework/math/Quadratures/Jacobi/quadrature_jacobi.h"

void
chi_math::QuadratureConical::Initialize_Conical_Product_Tri()
{
  QuadratureGaussLegendre legendre((QuadratureOrder)std::floor(((int)order_ + 1) / 2.0));
  QuadratureJacobi jacobiA(order_, /*alpha=*/1, /*beta=*/0);

  legendre.SetRange({0, 1});

  const size_t np = legendre.qpoints_.size();

  qpoints_.resize(np * np);
  weights_.resize(np * np);

  // Compute the conical product
  unsigned int gp = 0;
  for (unsigned int i = 0; i < np; i++)
    for (unsigned int j = 0; j < np; j++)
    {
      qpoints_[gp](0) = jacobiA.qpoints_[j](0);                                  // s[j];
      qpoints_[gp](1) = legendre.qpoints_[i](0) * (1. - jacobiA.qpoints_[j](0)); // r[i]*(1.-s[j]);
      weights_[gp] = legendre.weights_[i] * jacobiA.weights_[j];                 // A[i]*B[j];
      gp++;
    }
}

void
chi_math::QuadratureConical::Initialize_Conical_Product_Tet()
{
  QuadratureGaussLegendre legendre((QuadratureOrder)std::floor(((int)order_ + 1) / 2.0));
  QuadratureJacobi jacobiA(order_, /*alpha=*/1, /*beta=*/0);
  QuadratureJacobi jacobiB(order_, /*alpha=*/2, /*beta=*/0);

  legendre.SetRange({0, 1});

  const size_t np = legendre.qpoints_.size();

  qpoints_.resize(np * np * np);
  weights_.resize(np * np * np);

  double weight_sum = 0.0;
  unsigned int gp = 0;
  for (unsigned int i = 0; i < np; i++)
    for (unsigned int j = 0; j < np; j++)
      for (unsigned int k = 0; k < np; k++)
      {
        qpoints_[gp](0) = jacobiB.qpoints_[k](0);                                 // t[k];
        qpoints_[gp](1) = jacobiA.qpoints_[j](0) * (1. - jacobiB.qpoints_[k](0)); // s[j]*(1.-t[k]);
        qpoints_[gp](2) = legendre.qpoints_[i](0) * (1. - jacobiA.qpoints_[j](0)) *
                          (1. - jacobiB.qpoints_[k](0)); // r[i]*(1.-s[j])*(1.-t[k]);
        weights_[gp] =
          legendre.weights_[i] * jacobiA.weights_[j] * jacobiB.weights_[k]; // A[i]*B[j]*C[k];
        weight_sum += weights_[gp];
        gp++;
      }

  double w_scale = (1.0 / 6.0) / weight_sum;
  for (auto& v : weights_)
    v *= w_scale;
}
