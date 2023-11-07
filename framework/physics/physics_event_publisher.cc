#include "framework/physics/physics_event_publisher.h"
#include "framework/event_system/event.h"
#include "framework/event_system/system_wide_event_publisher.h"
#include "framework/event_system/event_codes.h"

#include "framework/physics/solver_base/solver.h"
#include "framework/physics/time_steppers/time_stepper.h"

namespace chi_physics
{

PhysicsEventPublisher::PhysicsEventPublisher() : chi::EventPublisher("Physics")
{
}

PhysicsEventPublisher&
PhysicsEventPublisher::GetInstance()
{
  static PhysicsEventPublisher singleton;
  return singleton;
}

void
PhysicsEventPublisher::PublishEvent(const chi::Event& event)
{
  chi::EventPublisher::PublishEvent(event);

  chi::SystemWideEventPublisher::GetInstance().PublishEvent(event);
}

void
PhysicsEventPublisher::SolverInitialize(Solver& solver)
{
  {
    const std::string event_name = "SolverPreInitialize";
    chi::ParameterBlock params;
    params.AddParameter("solver_name", solver.TextName());
    params.AddParameter("solver_handle", solver.StackID());

    PublishEvent(chi::Event(event_name, chi::GetStandardEventCode(event_name), params));
  }

  solver.Initialize();

  {
    const std::string event_name = "SolverInitialized";
    chi::ParameterBlock params;
    params.AddParameter("solver_name", solver.TextName());
    params.AddParameter("solver_handle", solver.StackID());
    params.AddParameter("time", solver.GetTimeStepper().Time());

    PublishEvent(chi::Event(event_name, chi::GetStandardEventCode(event_name), params));
  }
}

void
PhysicsEventPublisher::SolverExecute(Solver& solver)
{
  {
    const std::string event_name = "SolverPreExecution";
    chi::ParameterBlock params;
    params.AddParameter("solver_name", solver.TextName());
    params.AddParameter("solver_handle", solver.StackID());

    PublishEvent(chi::Event(event_name, chi::GetStandardEventCode(event_name), params));
  }

  solver.Execute();

  {
    const std::string event_name = "SolverExecuted";
    chi::ParameterBlock params;
    params.AddParameter("solver_name", solver.TextName());
    params.AddParameter("solver_handle", solver.StackID());

    PublishEvent(chi::Event(event_name, chi::GetStandardEventCode(event_name), params));
  }
}

void
PhysicsEventPublisher::SolverStep(Solver& solver)
{
  {
    const std::string event_name = "SolverPreStep";
    chi::ParameterBlock params;
    params.AddParameter("solver_name", solver.TextName());
    params.AddParameter("solver_handle", solver.StackID());

    PublishEvent(chi::Event(event_name, chi::GetStandardEventCode(event_name), params));
  }

  solver.Step();

  {
    const std::string event_name = "SolverStep";
    chi::ParameterBlock params;
    params.AddParameter("solver_name", solver.TextName());
    params.AddParameter("solver_handle", solver.StackID());

    PublishEvent(chi::Event(event_name, chi::GetStandardEventCode(event_name), params));
  }
}

void
PhysicsEventPublisher::SolverAdvance(Solver& solver)
{
  {
    const std::string event_name = "SolverPreAdvance";
    chi::ParameterBlock params;
    params.AddParameter("solver_name", solver.TextName());
    params.AddParameter("solver_handle", solver.StackID());

    PublishEvent(chi::Event(event_name, chi::GetStandardEventCode(event_name), params));
  }

  solver.Advance();

  {
    const std::string event_name = "SolverAdvanced";
    chi::ParameterBlock params;
    params.AddParameter("solver_name", solver.TextName());
    params.AddParameter("solver_handle", solver.StackID());
    params.AddParameter("time", solver.GetTimeStepper().Time());
    params.AddParameter("timestep_index", solver.GetTimeStepper().TimeStepIndex());

    PublishEvent(chi::Event(event_name, chi::GetStandardEventCode(event_name), params));
  }
}

} // namespace chi_physics
