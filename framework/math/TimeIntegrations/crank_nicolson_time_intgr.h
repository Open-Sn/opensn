#ifndef CHITECH_CRANK_NICOLSON_TIME_INTGR_H
#define CHITECH_CRANK_NICOLSON_TIME_INTGR_H

#include "theta_scheme_time_intgr.h"

namespace chi_math
{

class CrankNicolsonTimeIntegration : public ThetaSchemeTimeIntegration
{
public:
  static chi::InputParameters GetInputParameters();
  explicit CrankNicolsonTimeIntegration(
    const chi::InputParameters& params);
};

} // namespace chi_math

#endif // CHITECH_CRANK_NICOLSON_TIME_INTGR_H
