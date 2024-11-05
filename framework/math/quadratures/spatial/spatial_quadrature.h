// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/quadrature_order.h"
#include "framework/mesh/mesh.h"
#include "framework/object.h"
#include <vector>

namespace opensn
{

class SpatialQuadrature : public Object
{
protected:
  bool verbose_;
  QuadratureOrder order_;

  explicit SpatialQuadrature(const InputParameters& params)
    : Object(params),
      verbose_(params.ParamValue<bool>("verbose")),
      order_(static_cast<QuadratureOrder>(params.ParamValue<int>("order")))
  {
  }

  explicit SpatialQuadrature(QuadratureOrder order) : verbose_(false), order_(order) {}

public:
  std::vector<Vector3> qpoints;
  std::vector<double> weights;

  static InputParameters GetInputParameters();

  QuadratureOrder Order() { return order_; }
};

} // namespace opensn
