#include "opensn/framework/event_system/Event.h"

namespace chi
{

Event::Event(const std::string& name, int code) : name_(name), code_(code), params_()
{
}

Event::Event(const std::string& name, int code, const ParameterBlock& parameter_block)
  : name_(name), code_(code), params_(parameter_block)
{
}

const std::string&
Event::Name() const
{
  return name_;
}
int
Event::Code() const
{
  return code_;
}
const ParameterBlock&
Event::Parameters() const
{
  return params_;
}

} // namespace chi
