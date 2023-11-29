#include "framework/math/time_integrations/implicit_euler_time_intgr.h"

#include "framework/object_factory.h"

namespace opensn
{

OpenSnRegisterObject(math, ImplicitEulerTimeIntegration);

InputParameters
ImplicitEulerTimeIntegration::GetInputParameters()
{
  InputParameters params = ThetaSchemeTimeIntegration::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription("General implicit Euler Time Integration");
  params.SetDocGroup("DocTimeIntegrations");
  // clang-format on

  params.ChangeExistingParamToOptional("method", static_cast<int>(SteppingMethod::IMPLICIT_EULER));
  params.ChangeExistingParamToOptional("theta", 1.0);

  return params;
}

ImplicitEulerTimeIntegration::ImplicitEulerTimeIntegration(const InputParameters& params)
  : ThetaSchemeTimeIntegration(params)
{
}

} // namespace opensn
