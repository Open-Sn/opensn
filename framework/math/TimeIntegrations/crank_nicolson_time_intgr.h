#pragma once

#include "opensn/framework/math/TimeIntegrations/theta_scheme_time_intgr.h"

namespace chi_math
{

class CrankNicolsonTimeIntegration : public ThetaSchemeTimeIntegration
{
public:
  static chi::InputParameters GetInputParameters();
  explicit CrankNicolsonTimeIntegration(const chi::InputParameters& params);
};

} // namespace chi_math
