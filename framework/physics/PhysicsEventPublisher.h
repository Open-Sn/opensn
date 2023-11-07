#pragma once

#include "framework/event_system/EventPublisher.h"

namespace chi_physics
{

class Solver;

/**A singleton object that can be subscribed to for events.*/
class PhysicsEventPublisher : public chi::EventPublisher
{
public:
  static PhysicsEventPublisher& GetInstance();
  /// Deleted copy constructor
  PhysicsEventPublisher(const PhysicsEventPublisher&) = delete;
  PhysicsEventPublisher
  /// Deleted assignment operator
  operator=(const PhysicsEventPublisher&) = delete;

  void PublishEvent(const chi::Event& event) override;

  void SolverInitialize(Solver& solver);
  void SolverExecute(Solver& solver);
  void SolverStep(Solver& solver);
  void SolverAdvance(Solver& solver);

private:
  PhysicsEventPublisher();
};

} // namespace chi_physics
