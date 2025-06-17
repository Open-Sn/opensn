// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/event_system/physics_event_publisher.h"
#include "framework/event_system/event.h"
#include "framework/event_system/system_wide_event_publisher.h"
#include "physics/solvers/solver.h"
#include "physics/solvers/time_steppers/time_stepper.h"

namespace opensn
{

PhysicsEventPublisher::PhysicsEventPublisher() : EventPublisher("Physics")
{
}

PhysicsEventPublisher&
PhysicsEventPublisher::GetInstance()
{
  static PhysicsEventPublisher singleton;
  return singleton;
}

void
PhysicsEventPublisher::PublishEvent(const Event& event)
{
  EventPublisher::PublishEvent(event);

  SystemWideEventPublisher::GetInstance().PublishEvent(event);
}

void
PhysicsEventPublisher::SolverInitialize(Solver& solver)
{
  {
    const std::string event_name = "SolverPreInitialize";
    ParameterBlock params;
    params.AddParameter("solver_name", solver.GetName());
    params.AddParameter("solver_handle", solver.GetStackID());

    PublishEvent(Event(event_name, params));
  }

  solver.Initialize();

  {
    const std::string event_name = "SolverInitialized";
    ParameterBlock params;
    params.AddParameter("solver_name", solver.GetName());
    params.AddParameter("solver_handle", solver.GetStackID());
    params.AddParameter("time", solver.GetTimeStepper().GetTime());

    PublishEvent(Event(event_name, params));
  }
}

void
PhysicsEventPublisher::SolverExecute(Solver& solver)
{
  {
    const std::string event_name = "SolverPreExecution";
    ParameterBlock params;
    params.AddParameter("solver_name", solver.GetName());
    params.AddParameter("solver_handle", solver.GetStackID());

    PublishEvent(Event(event_name, params));
  }

  solver.Execute();

  {
    const std::string event_name = "SolverExecuted";
    ParameterBlock params;
    params.AddParameter("solver_name", solver.GetName());
    params.AddParameter("solver_handle", solver.GetStackID());

    PublishEvent(Event(event_name, params));
  }
}

void
PhysicsEventPublisher::SolverStep(Solver& solver)
{
  {
    const std::string event_name = "SolverPreStep";
    ParameterBlock params;
    params.AddParameter("solver_name", solver.GetName());
    params.AddParameter("solver_handle", solver.GetStackID());

    PublishEvent(Event(event_name, params));
  }

  solver.Step();

  {
    const std::string event_name = "SolverStep";
    ParameterBlock params;
    params.AddParameter("solver_name", solver.GetName());
    params.AddParameter("solver_handle", solver.GetStackID());

    PublishEvent(Event(event_name, params));
  }
}

void
PhysicsEventPublisher::SolverAdvance(Solver& solver)
{
  {
    const std::string event_name = "SolverPreAdvance";
    ParameterBlock params;
    params.AddParameter("solver_name", solver.GetName());
    params.AddParameter("solver_handle", solver.GetStackID());

    PublishEvent(Event(event_name, params));
  }

  solver.Advance();

  {
    const std::string event_name = "SolverAdvanced";
    ParameterBlock params;
    params.AddParameter("solver_name", solver.GetName());
    params.AddParameter("solver_handle", solver.GetStackID());
    params.AddParameter("time", solver.GetTimeStepper().GetTime());
    params.AddParameter("timestep_index", solver.GetTimeStepper().GetTimeStepIndex());

    PublishEvent(Event(event_name, params));
  }
}

} // namespace opensn
