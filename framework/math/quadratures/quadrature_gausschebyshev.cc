#include "framework/math/quadratures/quadrature_gausschebyshev.h"

#include "framework/object_factory.h"
#include <cmath>

#include "framework/runtime.h"
#include "framework/logging/log.h"

#define uint unsigned int
#define scint static_cast<int>

namespace opensn
{

RegisterChiObject(chi_math, QuadratureGaussChebyshev);

InputParameters
QuadratureGaussChebyshev::GetInputParameters()
{
  InputParameters params = Quadrature::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "Implementation of a Gauss-Chebyshev quadrature");
  params.SetDocGroup("LuaQuadrature");
  // clang-format on

  params.ChangeExistingParamToOptional("order", 0);

  params.AddOptionalParameter("N", 1, "Number of quadrature points.");

  return params;
}

QuadratureGaussChebyshev::QuadratureGaussChebyshev(const InputParameters& params)
  : Quadrature(params)
{
  const auto& assigned_params = params.ParametersAtAssignment();

  const int param_count = int(assigned_params.Has("order")) + int(assigned_params.Has("N"));
  ChiInvalidArgumentIf(param_count == 2, "Either \"order\" or \"N\" must be specified, not both");

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

QuadratureGaussChebyshev::QuadratureGaussChebyshev(unsigned int N, bool verbose)
  : Quadrature((QuadratureOrder)(2 * N - 1))
{
  Initialize(N);
}

void
QuadratureGaussChebyshev::Initialize(unsigned int N)
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
      Chi::log.Log() << "root[" << n << "]=" << qpoints_[n][0] << ", weight=" << weights_[n];
  }

  range_ = {-1, +1};
}

} // namespace opensn
