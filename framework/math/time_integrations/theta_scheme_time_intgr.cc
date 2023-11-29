#include "framework/math/time_integrations/theta_scheme_time_intgr.h"

#include "framework/object_factory.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(math, ThetaSchemeTimeIntegration);

InputParameters
ThetaSchemeTimeIntegration::GetInputParameters()
{
  InputParameters params = TimeIntegration::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription("Generalized theta-scheme");
  params.SetDocGroup("DocTimeIntegrations");
  // clang-format on

  params.ChangeExistingParamToOptional("method", static_cast<int>(SteppingMethod::THETA_SCHEME));

  params.AddRequiredParameter<double>("theta", "The theta parameter for a theta scheme");

  return params;
}

ThetaSchemeTimeIntegration::ThetaSchemeTimeIntegration(const InputParameters& params)
  : TimeIntegration(params), theta_(params.GetParamValue<double>("theta"))
{
}

double
ThetaSchemeTimeIntegration::ThetaFactor() const
{
  return theta_;
}

} // namespace opensn
