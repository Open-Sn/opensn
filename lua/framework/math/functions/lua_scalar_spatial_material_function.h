#pragma once

#include "framework/math/functions/scalar_spatial_material_function.h"

namespace opensnlua
{

class LuaScalarSpatialMaterialFunction : public opensn::ScalarSpatialMaterialFunction
{
public:
  static opensn::InputParameters GetInputParameters();
  explicit LuaScalarSpatialMaterialFunction(const opensn::InputParameters& params);
  double Evaluate(int mat_id, const opensn::Vector3& xyz) const override;

private:
  const std::string lua_function_name_;
};

} // namespace opensnlua
