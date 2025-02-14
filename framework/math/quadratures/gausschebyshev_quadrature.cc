// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/gausschebyshev_quadrature.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <cmath>

namespace opensn
{

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
