// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/solvers/time_steppers/constant_time_stepper.h"
#include "framework/object_factory.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(physics, ConstantTimeStepper);

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
