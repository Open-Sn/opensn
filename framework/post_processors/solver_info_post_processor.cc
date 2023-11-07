#include "framework/post_processors/solver_info_post_processor.h"

#include "framework/object_factory.h"

#include "framework/physics/solver_base/solver.h"
#include "framework/physics/time_steppers/time_stepper.h"
#include "framework/event_system/event.h"

#include <algorithm>

namespace chi
{

RegisterChiObject(chi, SolverInfoPostProcessor);

InputParameters
SolverInfoPostProcessor::GetInputParameters()
{
  InputParameters params = PostProcessor::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
"A post processor that can get basic info for any object based on "
"chi_physics::Solver. This solver's execution does not filter whether solver"
"events are for the relevant solver. This is done to avoid differing "
"time-histories.");
  // clang-format on
  params.SetDocGroup("doc_PostProcessors");

  params.AddRequiredParameter<size_t>("solver", "A handle to the solver.");
  params.AddRequiredParameterBlock(
    "info",
    "A parameter block, requiring at minimum the parameter \"name\" to pass to "
    "the solver. Note: each solver will have custom parameter blocks available "
    "to "
    "get a vast array of different pieces of information.");

  return params;
}

SolverInfoPostProcessor::SolverInfoPostProcessor(const InputParameters& params)
  : PostProcessor(params, PPType::SCALAR),
    solver_(Chi::GetStackItem<chi_physics::Solver>(
      Chi::object_stack, params.GetParamValue<size_t>("solver"), __FUNCTION__)),
    info_(params.GetParam("info"))
{
  const auto& param_assigned = params.ParametersAtAssignment();
  if (param_assigned.Has("solvername_filter"))
    solvername_filter_ = params.GetParamValue<std::string>("solvername_filter");
  else
    solvername_filter_ = solver_.TextName();
}

void
SolverInfoPostProcessor::Execute(const Event& event_context)
{
  value_ = solver_.GetInfoWithPreCheck(info_);
  SetType(FigureTypeFromValue(value_));

  const int event_code = event_context.Code();
  if (event_code == 32 /*SolverInitialized*/ or event_code == 38 /*SolverAdvanced*/)
  {
    TimeHistoryEntry entry{
      solver_.GetTimeStepper().TimeStepIndex(), solver_.GetTimeStepper().Time(), value_};
    time_history_.push_back(std::move(entry));
  }
}

} // namespace chi
