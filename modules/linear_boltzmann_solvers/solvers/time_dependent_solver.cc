// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/time_dependent_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/object_factory.h"
#include "framework/utils/hdf_utils.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <memory>
#include "framework/logging/log.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, TimeDependentSourceSolver);

InputParameters
TimeDependentSourceSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();

  params.SetGeneralDescription("Implementation of a time-dependent solver.");
  params.ChangeExistingParamToOptional("name", "TimeDependentSourceSolver");
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem", "An existing lbs problem");
  params.AddRequiredParameter<double>("dt", "Time step size.");
  params.ConstrainParameterRange("dt", AllowableRangeLowLimit::New(1.0e-18));
  params.AddRequiredParameter<double>("stop_time", "Time duration to run the solver.");
  params.ConstrainParameterRange("stop_time", AllowableRangeLowLimit::New(1.0e-16));
  params.AddRequiredParameter<double>("theta", "Time differencing scheme");
  params.ConstrainParameterRange("theta", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("theta", AllowableRangeHighLimit::New(1.0));
  return params;
}

std::shared_ptr<TimeDependentSourceSolver>
TimeDependentSourceSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<TimeDependentSourceSolver>("lbs::TimeDependentSourceSolver", params);
}

TimeDependentSourceSolver::TimeDependentSourceSolver(const InputParameters& params)
  : Solver(params),
    lbs_problem_(params.GetSharedPtrParam<Problem, LBSProblem>("problem")),
    dt_(params.GetParamValue<double>("dt")),
    theta_(params.GetParamValue<double>("theta")), 
    stop_time_(params.GetParamValue<double>("stop_time"))
{
  if (not lbs_problem_->IsTimeDependent())
    throw std::runtime_error(GetName() + ": Problem is not time dependent.");
}

void
TimeDependentSourceSolver::Initialize()
{
  CALI_CXX_MARK_SCOPE("TimeDependentSourceSolver::Initialize");

  lbs_problem_->Initialize();
}

void
TimeDependentSourceSolver::Execute()
{
  CALI_CXX_MARK_SCOPE("TimeDependentSourceSolver::Execute");

  unsigned int step = 0;
  double current_time = 0.0;
  auto& ags_solver = *lbs_problem_->GetAGSSolver();

  while (current_time < stop_time_)
  {
    const double target_time = std::min(current_time + dt_, stop_time_);
    const double step_dt = target_time - current_time;

    lbs_problem_->SetTimeStep(step_dt);
    lbs_problem_->SetTheta(theta_);

    // Simulation time is the new time level t_{n+1}
    lbs_problem_->SetSimulationTime(target_time);

    log.Log() << "\n*************** Time step #" << step
              << "  t = " << std::setprecision(6) << target_time
              << "  (from " << current_time
              << ", dt = " << step_dt
              << ", theta = " << theta_
              << ") ***************\n";

    ags_solver.Solve();
    lbs_problem_->UpdateFieldFunctions();
    lbs_problem_->UpdatePsiOld();

    // Advance current_time to the end of this step
    current_time = target_time;
    ++step;
  }

  // Ensure final simulation time matches requested stop time
  lbs_problem_->SetSimulationTime(stop_time_);
}


} // namespace opensn
