#ifndef CHITECH_IMPLICIT_EULER_TIME_INTGR_H
#define CHITECH_IMPLICIT_EULER_TIME_INTGR_H

#include "theta_scheme_time_intgr.h"

namespace chi_math
{

class ImplicitEulerTimeIntegration : public ThetaSchemeTimeIntegration
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ImplicitEulerTimeIntegration(
    const chi::InputParameters& params);
};

} // namespace chi_math

#endif // CHITECH_IMPLICIT_EULER_TIME_INTGR_H
