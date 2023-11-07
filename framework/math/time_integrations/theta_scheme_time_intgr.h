#pragma once

#include "framework/math/time_integrations/time_integration.h"

namespace chi_math
{

class ThetaSchemeTimeIntegration : public TimeIntegration
{
private:
  const double theta_;

public:
  static chi::InputParameters GetInputParameters();
  explicit ThetaSchemeTimeIntegration(const chi::InputParameters& params);

  double ThetaFactor() const;
};

} // namespace chi_math
