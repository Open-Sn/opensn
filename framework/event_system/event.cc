#include "framework/event_system/event.h"

namespace chi
{

Event::Event(const std::string& name, EventCode code) : name_(name), code_(code), params_()
{
}

Event::Event(const std::string& name, EventCode code, const ParameterBlock& parameter_block)
  : name_(name), code_(code), params_(parameter_block)
{
}

const std::string&
Event::Name() const
{
  return name_;
}

Event::EventCode
Event::Code() const
{
  return code_;
}

const ParameterBlock&
Event::Parameters() const
{
  return params_;
}

bool
Event::IsSolverEvent() const
{
  return (code_ >= SolverPreInitialize) && (code_ <= SolverAdvanced);
}

} // namespace chi
