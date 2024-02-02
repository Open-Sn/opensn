#include "framework/post_processors/post_processor.h"

#include "framework/physics/physics_event_publisher.h"
#include "framework/event_system/event_subscriber.h"
#include "framework/event_system/event.h"

#include "framework/logging/log.h"
#include "framework/object_factory.h"
#include <inttypes.h>

namespace opensn
{

OpenSnRegisterObjectParametersOnly(PostProcessor);

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
  const auto& user_params = params.ParametersAtAssignment();
  if (user_params.Has("initial_value"))
  {
    value_ = params.GetParam("initial_value");
    SetType(FigureTypeFromValue(value_));
  }
}

PPNumericFormat
PostProcessor::ConstructNumericFormat(const std::string& format_string)
{
  if (format_string == "fixed") return PPNumericFormat::FIXED;
  else if (format_string == "floating_point")
    return PPNumericFormat::FLOATING_POINT;
  else if (format_string == "scientific")
    return PPNumericFormat::SCIENTIFIC;
  else if (format_string == "general")
    return PPNumericFormat::GENERAL;
  else
    ChiLogicalError("Invalid numeric format string \"" + format_string + "\"");
}

const std::string&
PostProcessor::Name() const
{
  return name_;
}
PPType
PostProcessor::Type() const
{
  return type_;
}

PPNumericFormat
PostProcessor::NumericFormat() const
{
  return print_numeric_format_;
}

size_t
PostProcessor::NumericPrecision() const
{
  return print_precision_;
}

void
PostProcessor::PushOntoStack(std::shared_ptr<Object>& new_object)
{

  auto pp_ptr = std::dynamic_pointer_cast<PostProcessor>(new_object);
  ChiLogicalErrorIf(not pp_ptr, "Failure to cast new object to PostProcessor");

  postprocessor_stack.push_back(pp_ptr);
  new_object->SetStackID(postprocessor_stack.size() - 1);

  auto new_subscriber = std::dynamic_pointer_cast<EventSubscriber>(pp_ptr);

  ChiLogicalErrorIf(not new_subscriber, "Failure to cast PostProcessor to EventSubscriber");

  auto& publisher = PhysicsEventPublisher::GetInstance();
  publisher.AddSubscriber(new_subscriber);
}

void
PostProcessor::ReceiveEventUpdate(const Event& event)
{
  auto it = std::find(
    subscribed_events_for_execution_.begin(), subscribed_events_for_execution_.end(), event.Name());

  if (it != subscribed_events_for_execution_.end())
  {
    if (event.IsSolverEvent() and not solvername_filter_.empty())
    {
      if (event.Parameters().GetParamValue<std::string>("solver_name") != solvername_filter_)
        return;
    }

    Execute(event);
    if (log.GetVerbosity() >= 1)
      log.Log0Verbose1() << "Post processor \"" << Name()
                         << "\" executed on "
                            "event \""
                         << event.Name() << "\".";
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
  if (value.Type() == ParameterBlockType::BOOLEAN)
  {
    value_string = value.GetValue<bool>() ? "true" : "false";
  }
  else if (value.Type() == ParameterBlockType::FLOAT)
  {
    const auto dblval = value.GetValue<double>();
    char buffer[30];
    const auto numeric_format = NumericFormat();
    const size_t precision = NumericPrecision();
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
  else if (value.Type() == ParameterBlockType::STRING)
  {
    value_string = value.GetValue<std::string>();
  }
  else if (value.Type() == ParameterBlockType::INTEGER)
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
  if (type == PPType::SCALAR) return ConvertScalarValueToString(value);
  else if (type == PPType::VECTOR)
  {
    if (value.NumParameters() == 0) return "";
    const auto& first_entry = value.GetParam(0);
    const auto first_entry_type = first_entry.Type();

    ChiLogicalErrorIf(FigureTypeFromValue(first_entry) != PPType::SCALAR,
                      "The entries of the vector value of post-processor \"" + Name() +
                        "\" must all be SCALAR.");

    std::string output;
    for (const auto& entry : value)
    {
      ChiLogicalErrorIf(entry.Type() != first_entry_type,
                        "Mixed typed encountered in the vector values of post-processor \"" +
                          Name() + "\"");
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

  if (not value.HasValue() and value.NumParameters() == 0) return PPType::NO_VALUE;
  else if (IsScalar(value.Type()))
    return PPType::SCALAR;
  else if (value.Type() == ParameterBlockType::ARRAY)
  {
    if (value.NumParameters() == 0) return PPType::NO_VALUE;
    else
    {
      if (IsScalar(value.GetParam(0).Type())) return PPType::VECTOR;
      else
        return PPType::ARBITRARY;
    }
  }
  else if (value.Type() == ParameterBlockType::BLOCK)
    return PPType::ARBITRARY;
  else
    ChiLogicalError("Unsupported type");
}

void
PostProcessor::SetType(PPType type)
{
  type_ = type;
}

} // namespace opensn
