#include "time_integration.h"

#include "chi_log_exceptions.h"

#define scint static_cast<int>

namespace chi_math
{

chi::InputParameters TimeIntegration::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<int>("method",
                                   "Integer representing time stepping scheme");

  return params;
}

TimeIntegration::TimeIntegration(const chi::InputParameters& params)
  : ChiObject(params)
{
  const int method_option = params.GetParamValue<int>("method");
  if (method_option == scint(SteppingMethod::EXPLICIT_EULER))
    method_ = SteppingMethod::EXPLICIT_EULER;
  else if (method_option == scint(SteppingMethod::IMPLICIT_EULER))
    method_ = SteppingMethod::IMPLICIT_EULER;
  else if (method_option == scint(SteppingMethod::CRANK_NICOLSON))
    method_ = SteppingMethod::CRANK_NICOLSON;
  else if (method_option == scint(SteppingMethod::THETA_SCHEME))
    method_ = SteppingMethod::THETA_SCHEME;
  else
    ChiInvalidArgument("Unsupported Time Integration scheme");
}

SteppingMethod TimeIntegration::Method() const { return method_; }

} // namespace chi_math