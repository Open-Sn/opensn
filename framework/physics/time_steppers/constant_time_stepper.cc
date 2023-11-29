#include "framework/physics/time_steppers/constant_time_stepper.h"

#include "framework/object_factory.h"

namespace opensn
{

OpenSnRegisterObject(physics, ConstantTimeStepper);

InputParameters
ConstantTimeStepper::GetInputParameters()
{
  InputParameters params = TimeStepper::GetInputParameters();

  params.SetGeneralDescription("Timestep controller that does not dynamically change.");
  params.SetDocGroup("doc_TimeStepControllers");

  return params;
}

ConstantTimeStepper::ConstantTimeStepper(const InputParameters& params) : TimeStepper(params)
{
}

} // namespace opensn
