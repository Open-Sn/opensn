// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/gausschebyshev_quadrature.h"
#include "framework/logging/log.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include <cmath>

namespace opensn
{

OpenSnRegisterObjectInNamespace(squad, GaussChebyshevQuadrature);

InputParameters
GaussChebyshevQuadrature::GetInputParameters()
{
  InputParameters params = GaussQuadrature::GetInputParameters();

  params.SetGeneralDescription("Implementation of a Gauss-Chebyshev quadrature");

  params.SetDocGroup("LuaQuadrature");

  params.ChangeExistingParamToOptional("order", 0);
  params.ConstrainParameterRange("order", AllowableRangeLowHighLimit::New(0, 43));

  params.AddOptionalParameter("N", 1, "Number of quadrature points.");

  return params;
}

std::shared_ptr<GaussChebyshevQuadrature>
GaussChebyshevQuadrature::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<GaussChebyshevQuadrature>("squad::GaussChebyshevQuadrature", params);
}

GaussChebyshevQuadrature::GaussChebyshevQuadrature(const InputParameters& params)
  : GaussQuadrature(params)
{
  const int param_count = int(params.IsParameterValid("order")) + int(params.IsParameterValid("N"));
  OpenSnInvalidArgumentIf(param_count == 2,
                          "Either \"order\" or \"N\" must be specified, not both");

  if (params.IsParameterValid("order"))
  {
    const auto n = static_cast<unsigned int>(order_);
    Initialize(n);
  }
  else
  {
    const auto n = params.GetParamValue<unsigned int>("N");
    order_ = static_cast<QuadratureOrder>(std::min(n, 43u));
    Initialize(n);
  }
}

GaussChebyshevQuadrature::GaussChebyshevQuadrature(unsigned int N, bool verbose)
  : GaussQuadrature((QuadratureOrder)(2 * N - 1))
{
  Initialize(N);
}

void
GaussChebyshevQuadrature::Initialize(unsigned int N)
{
  if (verbose_)
    log.Log() << "Initializing Gauss-Chebyshev Quadrature with " << N << " q-points";

  const double pi_N = M_PI / N;
  for (unsigned int n = 0; n < N; ++n)
  {
    const double xn = -std::cos((2 * n + 1) * pi_N / 2.0);
    const double wn = pi_N;

    qpoints.emplace_back(xn);
    weights.emplace_back(wn);

    if (verbose_)
      log.Log() << "root[" << n << "]=" << qpoints[n][0] << ", weight=" << weights[n];
  }
  range_ = {-1, +1};
}

} // namespace opensn
