// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/time_integrations/theta_scheme_time_intgr.h"
#include "framework/object_factory.h"

namespace opensn
{

InputParameters
ThetaSchemeTimeIntegration::GetInputParameters()
{
  InputParameters params = TimeIntegration::GetInputParameters();

  params.SetGeneralDescription("Generalized theta-scheme");
  params.SetDocGroup("DocTimeIntegrations");

  params.ChangeExistingParamToOptional("method", static_cast<int>(SteppingMethod::THETA_SCHEME));

  params.AddRequiredParameter<double>("theta", "The theta parameter for a theta scheme");

  return params;
}

ThetaSchemeTimeIntegration::ThetaSchemeTimeIntegration(const InputParameters& params)
  : TimeIntegration(params), theta_(params.ParamValue<double>("theta"))
{
}

double
ThetaSchemeTimeIntegration::ThetaFactor() const
{
  return theta_;
}

} // namespace opensn
