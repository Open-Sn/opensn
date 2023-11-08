#include "framework/math/TimeIntegrations/theta_scheme_time_intgr.h"

#include "framework/ChiObjectFactory.h"

#define scint static_cast<int>

namespace chi_math
{

RegisterChiObject(chi_math, ThetaSchemeTimeIntegration);

chi::InputParameters
ThetaSchemeTimeIntegration::GetInputParameters()
{
  chi::InputParameters params = TimeIntegration::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription("Generalized theta-scheme");
  params.SetDocGroup("DocTimeIntegrations");
  // clang-format on

  params.ChangeExistingParamToOptional("method", scint(SteppingMethod::THETA_SCHEME));

  params.AddRequiredParameter<double>("theta", "The theta parameter for a theta scheme");

  return params;
}

ThetaSchemeTimeIntegration::ThetaSchemeTimeIntegration(const chi::InputParameters& params)
  : TimeIntegration(params), theta_(params.GetParamValue<double>("theta"))
{
}

double
ThetaSchemeTimeIntegration::ThetaFactor() const
{
  return theta_;
}

} // namespace chi_math
