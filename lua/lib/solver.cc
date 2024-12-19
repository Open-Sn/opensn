// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/solver.h"
#include "framework/event_system/physics_event_publisher.h"

using namespace opensn;

namespace opensnlua
{

void
SolverInitialize(std::shared_ptr<opensn::Solver> solver)
{
  PhysicsEventPublisher::GetInstance().SolverInitialize(*solver);
}

void
SolverExecute(std::shared_ptr<opensn::Solver> solver)
{
  PhysicsEventPublisher::GetInstance().SolverExecute(*solver);
}

void
SolverStep(std::shared_ptr<opensn::Solver> solver)
{
  PhysicsEventPublisher::GetInstance().SolverStep(*solver);
}

void
SolverAdvance(std::shared_ptr<opensn::Solver> solver)
{
  PhysicsEventPublisher::GetInstance().SolverAdvance(*solver);
}

} // namespace opensnlua
