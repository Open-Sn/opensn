// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/quadrature_order.h"
#include "framework/mesh/mesh.h"
#include "framework/parameters/input_parameters.h"
#include <vector>

namespace opensn
{

class SpatialQuadrature
{
protected:
  bool verbose_;
  QuadratureOrder order_;

  explicit SpatialQuadrature(const InputParameters& params)
    : verbose_(params.GetParamValue<bool>("verbose")),
      order_(static_cast<QuadratureOrder>(params.GetParamValue<int>("order")))
  {
  }

  explicit SpatialQuadrature(QuadratureOrder order) : verbose_(false), order_(order) {}

public:
  std::vector<Vector3> qpoints;
  std::vector<double> weights;

  static InputParameters GetInputParameters();

  QuadratureOrder GetOrder() { return order_; }
};

} // namespace opensn
