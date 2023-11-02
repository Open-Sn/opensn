#ifndef CHITECH_EVENT_H
#define CHITECH_EVENT_H

#include "parameters/parameter_block.h"

#include <string>

namespace chi
{

class Event
{
public:
  Event(const std::string& name, int code);
  Event(const std::string& name,
        int code,
        const ParameterBlock& parameter_block);
  const std::string& Name() const;
  int Code() const;
  const ParameterBlock& Parameters() const;

  virtual ~Event() = default;

protected:
  const std::string name_;
  const int code_ = 0;
  const ParameterBlock params_;
};

} // namespace chi

#endif // CHITECH_EVENT_H
