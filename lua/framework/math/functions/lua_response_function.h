#pragma once

#include "framework/math/functions/response_function.h"
#include "framework/mesh/mesh_vector.h"

namespace opensnlua
{

class LuaResponseFunction : public opensn::ResponseFunction
{
public:
  static opensn::InputParameters GetInputParameters();
  explicit LuaResponseFunction(const opensn::InputParameters& params);
  std::vector<double>
  Evaluate(int num_groups, const opensn::Vector3& xyz, int mat_id) const override;

private:
  const std::string lua_function_name_;
};

} // namespace opensnlua
