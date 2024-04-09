// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/time_steppers/time_stepper.h"

namespace opensn
{

/**Timestep controller that does not dynamically change.*/
class ConstantTimeStepper : public TimeStepper
{
public:
  static InputParameters GetInputParameters();
  explicit ConstantTimeStepper(const InputParameters& params);
};

} // namespace opensn
