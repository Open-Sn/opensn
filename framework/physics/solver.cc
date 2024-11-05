// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/physics/solver.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/physics/time_steppers/constant_time_stepper.h"
#include "framework/object_factory.h"

namespace opensn
{

InputParameters
Solver::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name",
    "A text name to associate with the solver. This name will be used "
    "in status messages and verbose iterative convergence monitors.");

  params.AddOptionalParameter("dt", 0.01, "Desired initial timestep size.");
  params.AddOptionalParameter("time", 0.0, "Current time of the solver.");
  params.AddOptionalParameter("start_time", 0.0, "Transient start-time if applicable.");
  params.AddOptionalParameter("end_time", 1.0, "Transient end-time if applicable.");
  params.AddOptionalParameter(
    "max_time_steps", -1, "Maximum number of timesteps to allow. Negative values disables this.");

  params.AddOptionalParameter("timestepper",
                              0,
                              "Handle to a timestepper. If not supplied then a ConstantTimeStepper "
                              "will be created.");

  params.ConstrainParameterRange("dt", AllowableRangeLowLimit::New(1.0e-12));

  return params;
}

Solver::Solver(std::string name)
  : timestepper_(InitTimeStepper(GetInputParameters())), name_(std::move(name))
{
}

Solver::Solver(std::string name, std::initializer_list<BasicOption> options)
  : basic_options_(options),
    timestepper_(InitTimeStepper(GetInputParameters())),
    name_(std::move(name))
{
}

Solver::Solver(const InputParameters& params)
  : Object(params),
    timestepper_(InitTimeStepper(params)),
    name_(params.ParamValue<std::string>("name"))
{
}

std::shared_ptr<TimeStepper>
Solver::InitTimeStepper(const InputParameters& params)
{
  const auto& user_params = params.ParametersAtAssignment();

  if (user_params.Has("timestepper"))
  {
    auto stepper = GetStackItemPtrAsType<TimeStepper>(
      object_stack, params.ParamValue<size_t>("timestepper"), __FUNCTION__);

    stepper->SetTimeStepSize(params.ParamValue<double>("dt"));
    stepper->SetTime(params.ParamValue<double>("time"));
    stepper->SetStartTime(params.ParamValue<double>("start_time"));
    stepper->SetEndTime(params.ParamValue<double>("end_time"));
    stepper->SetMaxTimeSteps(params.ParamValue<int>("max_time_steps"));

    return stepper;
  }
  else
  {
    auto& factory = ObjectFactory::Instance();

    const std::string obj_type = "physics::ConstantTimeStepper";
    auto valid_params = factory.RegisteredObjectParameters(obj_type);
    ParameterBlock custom_params;

    if (params.NumParameters() != 0)
    {
      custom_params.AddParameter(params.Param("dt"));
      custom_params.AddParameter(params.Param("time"));
      custom_params.AddParameter(params.Param("start_time"));
      custom_params.AddParameter(params.Param("end_time"));
      custom_params.AddParameter(params.Param("max_time_steps"));
    }

    valid_params.AssignParameters(custom_params);

    auto stepper = std::make_shared<ConstantTimeStepper>(valid_params);
    object_stack.push_back(stepper);
    stepper->SetStackID(object_stack.size() - 1);

    return stepper;
  }
}

std::string
Solver::Name() const
{
  return name_;
}

BasicOptions&
Solver::GetBasicOptions()
{
  return basic_options_;
}

const BasicOptions&
Solver::GetBasicOptions() const
{
  return basic_options_;
}

std::vector<std::shared_ptr<FieldFunctionGridBased>>&
Solver::FieldFunctions()
{
  return field_functions_;
}

TimeStepper&
Solver::GetTimeStepper()
{
  OpenSnLogicalErrorIf(not timestepper_, "Bad trouble: Timestepper not assigned.");
  return *timestepper_;
}

const TimeStepper&
Solver::GetTimeStepper() const
{
  OpenSnLogicalErrorIf(not timestepper_, "Bad trouble: Timestepper not assigned.");
  return *timestepper_;
}

const std::vector<std::shared_ptr<FieldFunctionGridBased>>&
Solver::FieldFunctions() const
{
  return field_functions_;
}

void
Solver::Initialize()
{
  log.Log() << "\"Initialize()\" method not defined for " << Name();
}

void
Solver::Execute()
{
  log.Log() << "\"Execute()\" method not defined for " << Name();
}

void
Solver::Step()
{
  log.Log() << "\"Step()\" method not defined for " << Name();
}

void
Solver::Advance()
{
  log.Log() << "\"Advance()\" method not defined for " << Name();
}

ParameterBlock
Solver::Info(const ParameterBlock& params) const
{
  return ParameterBlock{};
}

ParameterBlock
Solver::InfoWithPreCheck(const ParameterBlock& params) const
{
  if (not params.Has("name"))
  {
    log.LogAllWarning() << "Solver::GetInfo called without "
                           "\"name\" in the parameter list";
    return ParameterBlock{};
  }
  return Info(params);
}

void
Solver::SetProperties(const ParameterBlock& params)
{
  for (const auto& param : params)
  {
    const std::string param_name = param.Name();

    if (param_name == "dt")
      timestepper_->SetTimeStepSize(param.Value<double>());
    if (param_name == "time")
      timestepper_->SetTime(param.Value<double>());
    if (param_name == "start_time")
      timestepper_->SetStartTime(param.Value<double>());
    if (param_name == "end_time")
      timestepper_->SetEndTime(param.Value<double>());
    if (param_name == "max_time_steps")
      timestepper_->SetMaxTimeSteps(param.Value<int>());
    if (param_name == "dt_min")
      timestepper_->SetMinimumTimeStepSize(param.Value<int>());
  }
}

} // namespace opensn
