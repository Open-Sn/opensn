#include "framework/math/time_integrations/crank_nicolson_time_intgr.h"

#include "framework/object_factory.h"

#define scint static_cast<int>

namespace chi_math
{

RegisterChiObject(chi_math, CrankNicolsonTimeIntegration);

chi::InputParameters
CrankNicolsonTimeIntegration::GetInputParameters()
{
  chi::InputParameters params = ThetaSchemeTimeIntegration::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription("General Crank-Nicolson Time Integration");
  params.SetDocGroup("DocTimeIntegrations");
  // clang-format on

  params.ChangeExistingParamToOptional("method", scint(SteppingMethod::CRANK_NICOLSON));
  params.ChangeExistingParamToOptional("theta", 0.5);

  return params;
}

CrankNicolsonTimeIntegration::CrankNicolsonTimeIntegration(const chi::InputParameters& params)
  : ThetaSchemeTimeIntegration(params)
{
}

} // namespace chi_math
