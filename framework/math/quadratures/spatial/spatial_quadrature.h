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
      verbose_(params.GetParamValue<bool>("verbose")),
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
