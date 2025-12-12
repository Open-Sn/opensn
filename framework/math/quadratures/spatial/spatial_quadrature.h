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
public:
  QuadratureOrder GetOrder() { return order_; }

  std::vector<Vector3> qpoints;
  std::vector<double> weights;

protected:
  explicit SpatialQuadrature(const InputParameters& params)
    : verbose_(params.GetParamValue<bool>("verbose")),
      order_(static_cast<QuadratureOrder>(params.GetParamValue<int>("order")))
  {
  }

  explicit SpatialQuadrature(QuadratureOrder order) : verbose_(false), order_(order) {}

  bool verbose_;
  QuadratureOrder order_;

public:
  static InputParameters GetInputParameters();
};

} // namespace opensn
