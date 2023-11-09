#pragma once

#include "framework/parameters/parameter_block.h"

#include <string>

namespace chi
{

class Event
{
public:
  enum EventCode
  {
    Unknown = 0,
    ProgramStart = 1,
    ProgramExecuted = 2,
    SolverPreInitialize = 31,
    SolverInitialized = 32,
    SolverPreExecution = 33,
    SolverExecuted = 34,
    SolverPreStep = 35,
    SolverStep = 36,
    SolverPreAdvance = 37,
    SolverAdvanced = 38
  };

  Event(const std::string& name);
  Event(const std::string& name, const ParameterBlock& parameter_block);
  const std::string& Name() const;
  EventCode Code() const;
  const ParameterBlock& Parameters() const;

  virtual ~Event() = default;

  /**
   * Is this event a solver event
   * @return `true` is event is a solver event, `false` otherwise
   */
  bool IsSolverEvent() const;

protected:
  const std::string name_;
  const EventCode code_;
  const ParameterBlock params_;

public:
  /**
   * Gets the standard event code associated with the event name.
   *
   * \param event_name Event name
   * \return Event code (`ECode`) associated with the event name. If no code is found then
   * `Unknown`.
   */
  static Event::EventCode GetStandardCode(const std::string& event_name);
};

} // namespace chi
