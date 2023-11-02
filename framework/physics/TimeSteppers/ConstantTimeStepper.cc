#include "ConstantTimeStepper.h"

#include "ChiObjectFactory.h"

namespace chi_physics
{

RegisterChiObject(chi_physics, ConstantTimeStepper);

chi::InputParameters
ConstantTimeStepper::GetInputParameters()
{
  chi::InputParameters params = TimeStepper::GetInputParameters();

  params.SetGeneralDescription("Timestep controller that does not dynamically change.");
  params.SetDocGroup("doc_TimeStepControllers");

  return params;
}

ConstantTimeStepper::ConstantTimeStepper(const chi::InputParameters& params) : TimeStepper(params)
{
}

} // namespace chi_physics