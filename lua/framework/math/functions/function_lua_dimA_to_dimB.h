#pragma once

#include "framework/math/functions/function_dimA_to_dimB.h"

namespace opensnlua
{

class LuaDimAToDimB : public opensn::FunctionDimAToDimB
{
private:
  const std::string lua_function_name_;

public:
  static opensn::InputParameters GetInputParameters();

  explicit LuaDimAToDimB(const opensn::InputParameters& params);

  std::vector<double> Evaluate(const std::vector<double>& vals) const override;

  bool HasSlope() const override { return false; }
  bool HasCurvature() const override { return false; }
};

} // namespace opensnlua
