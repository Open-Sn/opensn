#pragma once

#include "framework/materials/material_property.h"
#include "framework/math/functions/function_dimA_to_dimB.h"

namespace opensn
{

/**General material base class for a scalar material property that
 * is possibly a function of position (x,y,z), time t, and maybe a set
 * of variables (v).*/
class MaterialPropertyScalarFuncXYZTV : public MaterialProperty
{
protected:
  const FunctionDimAToDimB& function_;
  const std::vector<std::string> dependent_variables_;

public:
  static InputParameters GetInputParameters();
  explicit MaterialPropertyScalarFuncXYZTV(const InputParameters& params);

  double Evaluate(const std::vector<double>& vars);
};

} // namespace opensn
