#pragma once

#include "theta_scheme_time_intgr.h"

namespace chi_math
{

class ImplicitEulerTimeIntegration : public ThetaSchemeTimeIntegration
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ImplicitEulerTimeIntegration(const chi::InputParameters& params);
};

} // namespace chi_math
