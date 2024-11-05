// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/quadrature_order.h"
#include "framework/mesh/mesh_vector.h"
#include "framework/object.h"
#include <vector>

namespace opensn
{

class GaussQuadrature : public Object
{
protected:
  std::pair<double, double> range_;
  bool verbose_;
  QuadratureOrder order_;

  explicit GaussQuadrature(const InputParameters& params)
    : Object(params),
      range_({0, 0}),
      verbose_(params.GetParamValue<bool>("verbose")),
      order_(static_cast<QuadratureOrder>(params.GetParamValue<int>("order")))
  {
  }

  explicit GaussQuadrature(QuadratureOrder order) : range_({0, 0}), verbose_(false), order_(order)
  {
  }

public:
  std::vector<Vector3> qpoints;
  std::vector<double> weights;

  static InputParameters GetInputParameters();

  QuadratureOrder Order() { return order_; }

  const std::pair<double, double>& Range() const { return range_; }

  void SetRange(const std::pair<double, double>& range);
};

} // namespace opensn
