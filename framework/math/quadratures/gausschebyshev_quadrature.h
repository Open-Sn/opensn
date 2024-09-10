// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/gauss_quadrature.h"
#include "framework/object_factory.h"

namespace opensn
{

class GaussChebyshevQuadrature : public GaussQuadrature
{
private:
  /**
   * Populates the abscissae and weights for a Gauss-Chebyshev quadrature given the number of
   * desired quadrature points.
   */
  void Initialize(unsigned int N);

public:
  explicit GaussChebyshevQuadrature(const InputParameters& params);

  /**
   * Populates the abscissae and weights for a Gauss-Chebyshev quadrature given the number of
   * desired quadrature points. The order of the quadrature will be 2N-1.
   */
  explicit GaussChebyshevQuadrature(unsigned int N, bool verbose = false);

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<GaussChebyshevQuadrature> Create(const ParameterBlock& params);
};

} // namespace opensn
