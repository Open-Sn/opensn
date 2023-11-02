#include "quadrature_gausschebyshev.h"

#include "ChiObjectFactory.h"
#include <cmath>

#include "chi_runtime.h"
#include "chi_log.h"

#define uint unsigned int
#define scint static_cast<int>

namespace chi_math
{

RegisterChiObject(chi_math, QuadratureGaussChebyshev);

chi::InputParameters QuadratureGaussChebyshev::GetInputParameters()
{
  chi::InputParameters params = Quadrature::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "Implementation of a Gauss-Chebyshev quadrature");
  params.SetDocGroup("LuaQuadrature");
  // clang-format on

  params.ChangeExistingParamToOptional("order", 0);

  params.AddOptionalParameter("N", 1, "Number of quadrature points.");

  return params;
}

QuadratureGaussChebyshev::QuadratureGaussChebyshev(
  const chi::InputParameters& params)
  : chi_math::Quadrature(params)
{
  const auto& assigned_params = params.ParametersAtAssignment();

  const int param_count =
    int(assigned_params.Has("order")) + int(assigned_params.Has("N"));
  ChiInvalidArgumentIf(param_count == 2,
                       "Either \"order\" or \"N\" must be specified, not both");

  if (assigned_params.Has("order"))
  {
    const uint N = static_cast<uint>(order_);
    Initialize(N);
  }
  else
  {
    const uint N = assigned_params.GetParamValue<uint>("N");
    order_ = static_cast<QuadratureOrder>(std::min(scint(N), 43));
    Initialize(N);
  }
}

// ###################################################################
/**Populates the abscissae and weights for a Gauss-Chebyshev
 * quadrature given the number of desired quadrature points. The
 * order of the quadrature will be 2N-1.*/
QuadratureGaussChebyshev::QuadratureGaussChebyshev(unsigned int N, bool verbose)
  : chi_math::Quadrature((QuadratureOrder)(2 * N - 1))
{
  Initialize(N);
}

// ###################################################################
/**Populates the abscissae and weights for a Gauss-Chebyshev
 * quadrature given the number of desired quadrature points.*/
void QuadratureGaussChebyshev::Initialize(unsigned int N)
{
  if (verbose_)
    Chi::log.Log() << "Initializing Gauss-Chebyshev Quadrature "
                      "with "
                   << N << " q-points";

  const double pi_N = M_PI / N;
  for (unsigned int n = 0; n < N; ++n)
  {
    const double xn = -std::cos((2 * n + 1) * pi_N / 2.0);
    const double wn = pi_N;

    qpoints_.emplace_back(xn);
    weights_.emplace_back(wn);

    if (verbose_)
      Chi::log.Log() << "root[" << n << "]=" << qpoints_[n][0]
                     << ", weight=" << weights_[n];
  }

  range_ = {-1, +1};
}

} // namespace chi_math