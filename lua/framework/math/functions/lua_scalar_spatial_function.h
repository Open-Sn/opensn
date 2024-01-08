#pragma once

#include "framework/math/functions/scalar_spatial_function.h"

namespace opensnlua
{

class LuaScalarSpatialFunction : public opensn::ScalarSpatialFunction
{
public:
  static opensn::InputParameters GetInputParameters();
  explicit LuaScalarSpatialFunction(const opensn::InputParameters& params);
  double Evaluate(const opensn::Vector3& xyz) const override;

private:
  const std::string lua_function_name_;
};

} // namespace opensnlua
