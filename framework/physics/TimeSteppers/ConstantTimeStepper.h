#pragma once

#include "TimeStepper.h"

namespace chi_physics
{

/**Timestep controller that does not dynamically change.*/
class ConstantTimeStepper : public TimeStepper
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ConstantTimeStepper(const chi::InputParameters& params);
};

} // namespace chi_physics


