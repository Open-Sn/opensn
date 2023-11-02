#ifndef CHI_MATH_QUADRATURE_GAUSS_CHEBYSHEV_H
#define CHI_MATH_QUADRATURE_GAUSS_CHEBYSHEV_H

#include "quadrature.h"

namespace chi_math
{

//######################################################### Class Def
/**Gauss-Chebyshev quadrature.*/
class QuadratureGaussChebyshev : public chi_math::Quadrature
{
public:
  static chi::InputParameters GetInputParameters();
  explicit QuadratureGaussChebyshev(const chi::InputParameters& params);

  explicit QuadratureGaussChebyshev(unsigned int N, bool verbose = false);

private:
  void Initialize(unsigned int N);
};

} // namespace chi_math
#endif // CHI_MATH_QUADRATURE_GAUSS_CHEBYSHEV_H