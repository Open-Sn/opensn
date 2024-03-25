#include "framework/materials/mat_prop_scalarfuncXYZTV.h"

#include "framework/object_factory.h"
#include "framework/runtime.h"

namespace opensn
{

InputParameters
MaterialPropertyScalarFuncXYZTV::GetInputParameters()
{
  InputParameters params = MaterialProperty::GetInputParameters();

  params.AddRequiredParameter<size_t>("function_handle",
                                      "Handle to a function to be used for "
                                      "evaluation of this material property.");
  params.AddOptionalParameterArray("dependent_variables",
                                   std::vector<std::string>{},
                                   "List of variable names denoting how the associated "
                                   "propery function will be called.");

  return params;
}

typedef FunctionDimAToDimB SFXYZV;

MaterialPropertyScalarFuncXYZTV::MaterialPropertyScalarFuncXYZTV(const InputParameters& params)
  : MaterialProperty(params),
    function_(GetStackItem<SFXYZV>(
      object_stack, params.GetParamValue<size_t>("function_handle"), __FUNCTION__)),
    dependent_variables_(params.GetParamVectorValue<std::string>("dependent_variables"))
{
  printf("Test eval: %g %g\n", 400.0, Evaluate({400.0}));
  for (const auto& dep_var : dependent_variables_)
    printf("%s\n", dep_var.c_str());
}

double
MaterialPropertyScalarFuncXYZTV::Evaluate(const std::vector<double>& vars)
{
  return function_.Evaluate(vars).front();
}

} // namespace opensn
