// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/spatial/spatial_quadrature.h"

namespace opensn
{

InputParameters
SpatialQuadrature::GetInputParameters()
{
  InputParameters params;

  params.SetGeneralDescription("Base class for spatial quadratures");

  params.AddRequiredParameter<int>("order", "Quadrature order.");

  params.AddOptionalParameter("verbose", false, "Enables verbose operations");

  return params;
}

} // namespace opensn
