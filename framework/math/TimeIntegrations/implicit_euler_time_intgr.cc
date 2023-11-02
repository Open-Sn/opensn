#include "implicit_euler_time_intgr.h"

#include "ChiObjectFactory.h"

#define scint static_cast<int>

namespace chi_math
{

RegisterChiObject(chi_math, ImplicitEulerTimeIntegration);

chi::InputParameters ImplicitEulerTimeIntegration::GetInputParameters()
{
  chi::InputParameters params =
    ThetaSchemeTimeIntegration::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription("General implicit Euler Time Integration");
  params.SetDocGroup("DocTimeIntegrations");
  // clang-format on

  params.ChangeExistingParamToOptional("method",
                                       scint(SteppingMethod::IMPLICIT_EULER));
  params.ChangeExistingParamToOptional("theta", 1.0);

  return params;
}

ImplicitEulerTimeIntegration::ImplicitEulerTimeIntegration(
  const chi::InputParameters& params)
  : ThetaSchemeTimeIntegration(params)
{
}

} // namespace chi_math