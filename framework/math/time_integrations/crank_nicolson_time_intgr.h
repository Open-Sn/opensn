#pragma once

#include "framework/math/time_integrations/theta_scheme_time_intgr.h"

namespace opensn
{

class CrankNicolsonTimeIntegration : public ThetaSchemeTimeIntegration
{
public:
  static InputParameters GetInputParameters();
  explicit CrankNicolsonTimeIntegration(const InputParameters& params);
};

} // namespace opensn
