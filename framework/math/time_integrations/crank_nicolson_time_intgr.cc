// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/time_integrations/crank_nicolson_time_intgr.h"
#include "framework/object_factory.h"

namespace opensn
{

InputParameters
CrankNicolsonTimeIntegration::GetInputParameters()
{
  InputParameters params = ThetaSchemeTimeIntegration::GetInputParameters();

  params.SetGeneralDescription("General Crank-Nicolson Time Integration");
  params.SetDocGroup("DocTimeIntegrations");

  params.ChangeExistingParamToOptional("method", static_cast<int>(SteppingMethod::CRANK_NICOLSON));
  params.ChangeExistingParamToOptional("theta", 0.5);

  return params;
}

CrankNicolsonTimeIntegration::CrankNicolsonTimeIntegration(const InputParameters& params)
  : ThetaSchemeTimeIntegration(params)
{
}

} // namespace opensn
