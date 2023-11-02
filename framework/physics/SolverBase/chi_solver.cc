#include "chi_solver.h"

#include "utils/chi_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "physics/TimeSteppers/ConstantTimeStepper.h"

#include "ChiObjectFactory.h"

namespace chi_physics
{

/**Returns the input parameters.*/
chi::InputParameters
Solver::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

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

  using namespace chi_data_types;
  params.ConstrainParameterRange("dt", AllowableRangeLowLimit::New(1.0e-12));

  return params;
}

Solver::Solver(std::string in_text_name)
  : timestepper_(InitTimeStepper(GetInputParameters())), text_name_(std::move(in_text_name))
{
}

Solver::Solver(std::string in_text_name, std::initializer_list<BasicOption> in_options)
  : basic_options_(in_options),
    timestepper_(InitTimeStepper(GetInputParameters())),
    text_name_(std::move(in_text_name))
{
}

Solver::Solver(const chi::InputParameters& params)
  : ChiObject(params),
    timestepper_(InitTimeStepper(params)),
    text_name_(params.GetParamValue<std::string>("name"))
{
}

std::shared_ptr<TimeStepper>
Solver::InitTimeStepper(const chi::InputParameters& params)
{
  const auto& user_params = params.ParametersAtAssignment();

  if (user_params.Has("timestepper"))
  {
    auto stepper = Chi::GetStackItemPtrAsType<TimeStepper>(
      /*stack=*/Chi::object_stack,
      /*handle=*/params.GetParamValue<size_t>("timestepper"),
      /*calling_function_name=*/__FUNCTION__);

    stepper->SetTimeStepSize(params.GetParamValue<double>("dt"));
    stepper->SetTime(params.GetParamValue<double>("time"));
    stepper->SetStartTime(params.GetParamValue<double>("start_time"));
    stepper->SetEndTime(params.GetParamValue<double>("end_time"));
    stepper->SetMaxTimeSteps(params.GetParamValue<int>("max_time_steps"));

    return stepper;
  }
  else
  {
    auto& factory = ChiObjectFactory::GetInstance();

    const std::string obj_type = "chi_physics::ConstantTimeStepper";
    auto valid_params = factory.GetRegisteredObjectParameters(obj_type);
    chi::ParameterBlock custom_params;

    if (params.NumParameters() != 0)
    {
      custom_params.AddParameter(params.GetParam("dt"));
      custom_params.AddParameter(params.GetParam("time"));
      custom_params.AddParameter(params.GetParam("start_time"));
      custom_params.AddParameter(params.GetParam("end_time"));
      custom_params.AddParameter(params.GetParam("max_time_steps"));
    }

    valid_params.AssignParameters(custom_params);

    auto stepper = std::make_shared<ConstantTimeStepper>(valid_params);
    Chi::object_stack.push_back(stepper);
    stepper->SetStackID(Chi::object_stack.size() - 1);

    return stepper;
  }
}

std::string
Solver::TextName() const
{
  return text_name_;
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
Solver::GetFieldFunctions()
{
  return field_functions_;
}

TimeStepper&
Solver::GetTimeStepper()
{
  ChiLogicalErrorIf(not timestepper_, "Bad trouble: Timestepper not assigned.");
  return *timestepper_;
}

const TimeStepper&
Solver::GetTimeStepper() const
{
  ChiLogicalErrorIf(not timestepper_, "Bad trouble: Timestepper not assigned.");
  return *timestepper_;
}

const std::vector<std::shared_ptr<FieldFunctionGridBased>>&
Solver::GetFieldFunctions() const
{
  return field_functions_;
}

void
Solver::Initialize()
{
  Chi::log.Log() << "\"Initialize()\" method not defined for " << TextName();
}

void
Solver::Execute()
{
  Chi::log.Log() << "\"Execute()\" method not defined for " << TextName();
}

void
Solver::Step()
{
  Chi::log.Log() << "\"Step()\" method not defined for " << TextName();
}

void
Solver::Advance()
{
  Chi::log.Log() << "\"Advance()\" method not defined for " << TextName();
}

chi::ParameterBlock
Solver::GetInfo(const chi::ParameterBlock& params) const
{
  return chi::ParameterBlock{};
}

chi::ParameterBlock
Solver::GetInfoWithPreCheck(const chi::ParameterBlock& params) const
{
  if (not params.Has("name"))
  {
    Chi::log.LogAllWarning() << "chi_physics::Solver::GetInfo called without "
                                "\"name\" in the parameter list";
    return chi::ParameterBlock{};
  }
  return GetInfo(params);
}

/**\addtogroup SolverBase
 *
 * \section Properties Properties that can be set
 * The following properties can be set via the lua call
 * `chi_lua::chiSolverSetProperties`
 * \copydoc chi_physics::Solver::SetProperties*/

/**
Base solver settable properties:
* - `dt`, Timestep size
* - `time`, Current time
* */
void
Solver::SetProperties(const chi::ParameterBlock& params)
{
  for (const auto& param : params)
  {
    const std::string param_name = param.Name();

    if (param_name == "dt") timestepper_->SetTimeStepSize(param.GetValue<double>());
    if (param_name == "time") timestepper_->SetTime(param.GetValue<double>());
    if (param_name == "start_time") timestepper_->SetStartTime(param.GetValue<double>());
    if (param_name == "end_time") timestepper_->SetEndTime(param.GetValue<double>());
    if (param_name == "max_time_steps") timestepper_->SetMaxTimeSteps(param.GetValue<int>());
    if (param_name == "dt_min") timestepper_->SetMinimumTimeStepSize(param.GetValue<int>());
  }
}

} // namespace chi_physics
