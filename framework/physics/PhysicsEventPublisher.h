#ifndef CHITECH_PHYSICSEVENTPUBLISHER_H
#define CHITECH_PHYSICSEVENTPUBLISHER_H

#include "event_system/EventPublisher.h"

namespace chi_physics
{

class Solver;

/**A singleton object that can be subscribed to for events.*/
class PhysicsEventPublisher : public chi::EventPublisher
{
public:
  static PhysicsEventPublisher& GetInstance();
  PhysicsEventPublisher(const PhysicsEventPublisher&) = delete; // Deleted copy constructor
  PhysicsEventPublisher
  operator=(const PhysicsEventPublisher&) = delete; // Deleted assignment operator

  void PublishEvent(const chi::Event& event) override;

  void SolverInitialize(Solver& solver);
  void SolverExecute(Solver& solver);
  void SolverStep(Solver& solver);
  void SolverAdvance(Solver& solver);

private:
  PhysicsEventPublisher();
};

} // namespace chi_physics

#endif // CHITECH_PHYSICSEVENTPUBLISHER_H
