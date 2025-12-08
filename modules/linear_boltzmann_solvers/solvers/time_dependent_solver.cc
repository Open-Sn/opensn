// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/time_dependent_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/object_factory.h"
#include "framework/utils/hdf_utils.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <stdexcept>
#include <memory>

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
  params.AddOptionalParameter<double>("dt", 0.1, "Time step size.");
  params.ConstrainParameterRange("dt", AllowableRangeLowLimit::New(1.0e-18));
  params.AddOptionalParameter<double>("stop_time", 1.0, "Time duration to run the solver.");
  params.ConstrainParameterRange("stop_time", AllowableRangeLowLimit::New(1.0e-16));
  params.AddOptionalParameter<double>("theta", 1.0, "Time differencing scheme");
  params.ConstrainParameterRange("theta", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("theta", AllowableRangeHighLimit::New(1.0));
  params.AddOptionalParameter("verbose", true, "Verbose logging");
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
    stop_time_(params.GetParamValue<double>("stop_time")),
    verbose_(params.GetParamValue<bool>("verbose"))
{
  if (not lbs_problem_->IsTimeDependent())
    throw std::runtime_error(GetName() + ": Problem is not time dependent.");
}

void
TimeDependentSourceSolver::Initialize()
{
  CALI_CXX_MARK_SCOPE("TimeDependentSourceSolver::Initialize");

  lbs_problem_->Initialize();
  ags_solver_ = lbs_problem_->GetAGSSolver();
}

void
TimeDependentSourceSolver::Execute()
{
  CALI_CXX_MARK_SCOPE("TimeDependentSourceSolver::Execute");

  while (current_time_ < stop_time_)
  {
    const double target_time = std::min(current_time_ + dt_, stop_time_);
    const double step_dt = target_time - current_time_;
    SetTimeStep(step_dt);
    Advance();
  }
}

void
TimeDependentSourceSolver::Advance()
{
  const double t_start = current_time_;
  const double t_end = current_time_ + dt_;

  lbs_problem_->SetTimeStep(dt_);
  lbs_problem_->SetTheta(theta_);
  lbs_problem_->SetTime(t_start);

  if (verbose_)
  {
    log.Log() << "\n*************** Time step #" << step_ << "  t = " << std::setprecision(6)
              << t_end << "  (from = " << t_start << ", dt = " << dt_ << ", theta = " << theta_
              << ") ***************\n";
  }

  ags_solver_->Solve();
  lbs_problem_->UpdateFieldFunctions();
  lbs_problem_->UpdatePsiOld();

  current_time_ = t_end;
  lbs_problem_->SetTime(current_time_);
  ++step_;
}

void
TimeDependentSourceSolver::SetTimeStep(double dt)
{
  if (dt < 0.0)
    throw std::runtime_error(GetName() + " dt must be non-negative");
  dt_ = dt;
}

void
TimeDependentSourceSolver::SetTheta(double theta)
{
  if (theta < 0.0 or theta > 1.0)
    throw std::runtime_error(GetName() + " theta must be between 0.0 and 1.0.");
  theta_ = theta;
}

} // namespace opensn
