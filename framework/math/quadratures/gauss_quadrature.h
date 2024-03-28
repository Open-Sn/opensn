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
  QuadratureOrder order_;
  std::vector<Vector3> qpoints_;
  std::vector<double> weights_;

  static InputParameters GetInputParameters();

  const std::pair<double, double>& GetRange() const { return range_; }

  void SetRange(const std::pair<double, double>& range);
};

} // namespace opensn
