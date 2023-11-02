#ifndef CHITECH_CONSTANTTIMESTEPPER_H
#define CHITECH_CONSTANTTIMESTEPPER_H

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

}

#endif // CHITECH_CONSTANTTIMESTEPPER_H
