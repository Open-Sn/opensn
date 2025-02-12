// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/post_processors/post_processor.h"
#include "framework/event_system/physics_event_publisher.h"
#include "framework/event_system/event_subscriber.h"
#include "framework/event_system/event.h"
#include "framework/logging/log.h"
#include <cinttypes>

namespace opensn
{

InputParameters
PostProcessor::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("Base class for Post-Processors. For more general"
                               "information see \\ref doc_PostProcessors");
  params.SetDocGroup("doc_PostProcessors");

  params.AddRequiredParameter<std::string>(
    "name",
    "Name of the post processor. This name will be used in many places so make "
    "sure it's a useful name.");
  params.AddOptionalParameterArray(
    "execute_on",
    std::vector<std::string>{
      "SolverInitialized", "SolverAdvanced", "SolverExecuted", "ProgramExecuted"},
    "List of events at which the post-processor will execute.");

  params.AddOptionalParameterArray(
    "print_on",
    std::vector<std::string>{
      "SolverInitialized", "SolverAdvanced", "SolverExecuted", "ProgramExecuted"},
    "List of events at which the post-processor will print. Make sure that "
    "these "
    "events are also set for the `PostProcessorPrinter` otherwise it wont "
    "print.");

  params.AddOptionalParameterBlock("initial_value", ParameterBlock{}, "An initial value.");

  params.AddOptionalParameter("print_numeric_format", "general", "Numeric format to use.");

  params.ConstrainParameterRange(
    "print_numeric_format",
    AllowableRangeList::New({"fixed", "floating_point", "scientific", "general"}));

  params.AddOptionalParameter(
    "print_precision", 6, "Number of digits to display after decimal point");

  params.AddOptionalParameter("solvername_filter",
                              "",
                              "Controls update events to only execute on the relevant solver's"
                              "event calls.");

  return params;
}

PostProcessor::PostProcessor(const InputParameters& params, PPType type)
  : Object(params),
    name_(params.GetParamValue<std::string>("name")),
    subscribed_events_for_execution_(params.GetParamVectorValue<std::string>("execute_on")),
    subscribed_events_for_printing_(params.GetParamVectorValue<std::string>("print_on")),
    type_(type),
    print_numeric_format_(
      ConstructNumericFormat(params.GetParamValue<std::string>("print_numeric_format"))),
    print_precision_(params.GetParamValue<size_t>("print_precision")),
    solvername_filter_(params.GetParamValue<std::string>("solvername_filter"))
{
  if (params.IsParameterValid("initial_value"))
  {
    value_ = params.GetParam("initial_value");
    SetType(FigureTypeFromValue(value_));
  }
}

PPNumericFormat
PostProcessor::ConstructNumericFormat(const std::string& format_string)
{
  if (format_string == "fixed")
    return PPNumericFormat::FIXED;
  else if (format_string == "floating_point")
    return PPNumericFormat::FLOATING_POINT;
  else if (format_string == "scientific")
    return PPNumericFormat::SCIENTIFIC;
  else if (format_string == "general")
    return PPNumericFormat::GENERAL;
  else
    OpenSnLogicalError("Invalid numeric format string \"" + format_string + "\"");
}

const std::string&
PostProcessor::GetName() const
{
  return name_;
}
PPType
PostProcessor::GetType() const
{
  return type_;
}

PPNumericFormat
PostProcessor::GetNumericFormat() const
{
  return print_numeric_format_;
}

size_t
PostProcessor::GetNumericPrecision() const
{
  return print_precision_;
}

void
PostProcessor::PushOntoStack(std::shared_ptr<Object> new_object)
{

  auto pp_ptr = std::dynamic_pointer_cast<PostProcessor>(new_object);
  OpenSnLogicalErrorIf(not pp_ptr, "Failure to cast new object to PostProcessor");

  postprocessor_stack.push_back(pp_ptr);
  new_object->SetStackID(postprocessor_stack.size() - 1);

  auto new_subscriber = std::dynamic_pointer_cast<EventSubscriber>(pp_ptr);

  OpenSnLogicalErrorIf(not new_subscriber, "Failure to cast PostProcessor to EventSubscriber");

  auto& publisher = PhysicsEventPublisher::GetInstance();
  publisher.AddSubscriber(new_subscriber);
}

void
PostProcessor::ReceiveEventUpdate(const Event& event)
{
  auto it = std::find(subscribed_events_for_execution_.begin(),
                      subscribed_events_for_execution_.end(),
                      event.GetName());

  if (it != subscribed_events_for_execution_.end())
  {
    if (event.IsSolverEvent() and not solvername_filter_.empty())
    {
      if (event.Parameters().GetParamValue<std::string>("solver_name") != solvername_filter_)
        return;
    }

    Execute(event);
    if (log.GetVerbosity() >= 1)
      log.Log0Verbose1() << "Post processor \"" << GetName()
                         << "\" executed on "
                            "event \""
                         << event.GetName() << "\".";
  }
}

const ParameterBlock&
PostProcessor::GetValue() const
{
  return value_;
}

const std::vector<PostProcessor::TimeHistoryEntry>&
PostProcessor::GetTimeHistory() const
{
  return time_history_;
}

const std::vector<std::string>&
PostProcessor::PrintScope() const
{
  return subscribed_events_for_printing_;
}

std::string
PostProcessor::ConvertScalarValueToString(const ParameterBlock& value) const
{
  std::string value_string;
  if (value.GetType() == ParameterBlockType::BOOLEAN)
  {
    value_string = value.GetValue<bool>() ? "true" : "false";
  }
  else if (value.GetType() == ParameterBlockType::FLOAT)
  {
    const auto dblval = value.GetValue<double>();
    char buffer[30];
    const auto numeric_format = GetNumericFormat();
    const size_t precision = GetNumericPrecision();
    if (numeric_format == PPNumericFormat::SCIENTIFIC)
    {
      const std::string format_spec = "%." + std::to_string(precision) + "e";
      snprintf(buffer, 30, format_spec.c_str(), dblval);
    }
    else if (numeric_format == PPNumericFormat::FLOATING_POINT)
    {
      const std::string format_spec = "%." + std::to_string(precision) + "f";
      snprintf(buffer, 30, format_spec.c_str(), dblval);
    }
    else // GENERAL
    {
      if (dblval < 1.0e-4)
      {
        const std::string format_spec = "%." + std::to_string(precision) + "e";
        snprintf(buffer, 30, format_spec.c_str(), dblval);
      }
      else if (dblval >= 1.0e-4 and dblval < 1.0e6)
      {
        const std::string format_spec = "%." + std::to_string(precision) + "f";
        snprintf(buffer, 30, format_spec.c_str(), dblval);
      }
      else
      {
        const std::string format_spec = "%." + std::to_string(precision) + "e";
        snprintf(buffer, 30, format_spec.c_str(), dblval);
      }
    } // if num_format

    value_string = buffer;
  }
  else if (value.GetType() == ParameterBlockType::STRING)
  {
    value_string = value.GetValue<std::string>();
  }
  else if (value.GetType() == ParameterBlockType::INTEGER)
  {
    const auto intval = value.GetValue<int64_t>();
    char buffer[30];
    snprintf(buffer, 30, "%" PRId64, intval);
    value_string = buffer;
  }

  return value_string;
}

std::string
PostProcessor::ConvertValueToString(const ParameterBlock& value) const
{
  const PPType type = FigureTypeFromValue(value);
  if (type == PPType::SCALAR)
    return ConvertScalarValueToString(value);
  else if (type == PPType::VECTOR)
  {
    if (value.GetNumParameters() == 0)
      return "";
    const auto& first_entry = value.GetParam(0);
    const auto first_entry_type = first_entry.GetType();

    OpenSnLogicalErrorIf(FigureTypeFromValue(first_entry) != PPType::SCALAR,
                         "The entries of the vector value of post-processor \"" + GetName() +
                           "\" must all be SCALAR.");

    std::string output;
    for (const auto& entry : value)
    {
      OpenSnLogicalErrorIf(entry.GetType() != first_entry_type,
                           "Mixed typed encountered in the vector values of post-processor \"" +
                             GetName() + "\"");
      output.append(ConvertScalarValueToString(entry) + " ");
    }

    return output;
  }
  else
  {
    std::string outstr;
    value.RecursiveDumpToString(outstr);
    std::replace(outstr.begin(), outstr.end(), '\n', ' ');
    return outstr;
  }
}

PPType
PostProcessor::FigureTypeFromValue(const ParameterBlock& value)
{
  const std::vector<ParameterBlockType> scalar_types = {ParameterBlockType::BOOLEAN,
                                                        ParameterBlockType::FLOAT,
                                                        ParameterBlockType::STRING,
                                                        ParameterBlockType::INTEGER};

  /**Lambda to check if this is a scalar*/
  auto IsScalar = [&scalar_types](const ParameterBlockType& block_type)
  { return std::find(scalar_types.begin(), scalar_types.end(), block_type) != scalar_types.end(); };

  if (not value.HasValue() and value.GetNumParameters() == 0)
    return PPType::NO_VALUE;
  else if (IsScalar(value.GetType()))
    return PPType::SCALAR;
  else if (value.GetType() == ParameterBlockType::ARRAY)
  {
    if (value.GetNumParameters() == 0)
      return PPType::NO_VALUE;
    else
    {
      if (IsScalar(value.GetParam(0).GetType()))
        return PPType::VECTOR;
      else
        return PPType::ARBITRARY;
    }
  }
  else if (value.GetType() == ParameterBlockType::BLOCK)
    return PPType::ARBITRARY;
  else
    OpenSnLogicalError("Unsupported type");
}

void
PostProcessor::SetType(PPType type)
{
  type_ = type;
}

} // namespace opensn
