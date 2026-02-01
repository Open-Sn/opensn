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
    stop_time_(params.GetParamValue<double>("stop_time")),
    verbose_(params.GetParamValue<bool>("verbose"))
{
  if (auto do_problem = std::dynamic_pointer_cast<DiscreteOrdinatesProblem>(lbs_problem_))
  {
    do_problem->EnableTimeDependentMode();
  }

  lbs_problem_->SetTimeStep(params.GetParamValue<double>("dt"));
  lbs_problem_->SetTheta(params.GetParamValue<double>("theta"));
}

void
TimeDependentSourceSolver::Initialize()
{
  CALI_CXX_MARK_SCOPE("TimeDependentSourceSolver::Initialize");

  lbs_problem_->Initialize();
  ags_solver_ = lbs_problem_->GetAGSSolver();
  initialized_ = true;
}

void
TimeDependentSourceSolver::Execute()
{
  if (not initialized_)
    throw std::runtime_error(GetName() + ": Initialize must be called before Execute.");
  const double t0 = current_time_;
  const double tf = stop_time_;

  if (tf < t0)
    throw std::runtime_error(GetName() + ": stop_time must be >= current_time");

  // Tolerance to account for floating-point drift from repeated additions/subtractions.
  // Use a small multiple of machine epsilon scaled to the time magnitude.
  const double tol =
    64.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(tf), std::abs(t0)});

  while (true)
  {
    const double remaining = tf - current_time_;

    if (remaining <= tol)
    {
      current_time_ = tf;
      break;
    }

    if (pre_advance_callback_)
      pre_advance_callback_();

    const double dt = lbs_problem_->GetTimeStep();
    if (dt <= 0.0)
      throw std::runtime_error(GetName() + ": dt must be positive");
    const double step_dt = (remaining < dt) ? remaining : dt;

    lbs_problem_->SetTimeStep(step_dt);
    lbs_problem_->SetTime(current_time_);

    Advance();

    if (post_advance_callback_)
      post_advance_callback_();

    if (std::abs(tf - current_time_) <= tol)
      current_time_ = tf;
  }
}

void
TimeDependentSourceSolver::Advance()
{
  if (not initialized_)
    throw std::runtime_error(GetName() + ": Initialize must be called before Advance.");
  const double dt = lbs_problem_->GetTimeStep();
  const double theta = lbs_problem_->GetTheta();

  if (verbose_)
  {
    log.Log() << "\n*************** Time step #" << step_ << "  t = " << std::setprecision(6)
              << current_time_ + dt << "  (from = " << current_time_ << ", dt = " << dt
              << ", theta = " << theta << ") ***************\n";
  }

  ags_solver_->Solve();
  lbs_problem_->UpdateFieldFunctions();
  if (IsBalanceEnabled())
    if (auto do_problem = std::dynamic_pointer_cast<DiscreteOrdinatesProblem>(lbs_problem_))
      ComputeBalance(*do_problem);
  lbs_problem_->UpdatePsiOld();

  current_time_ += dt;
  lbs_problem_->SetTime(current_time_);
  ++step_;
}

void
TimeDependentSourceSolver::SetTimeStep(double dt)
{
  lbs_problem_->SetTimeStep(dt);
}

void
TimeDependentSourceSolver::SetTheta(double theta)
{
  lbs_problem_->SetTheta(theta);
}

void
TimeDependentSourceSolver::SetPreAdvanceCallback(std::function<void()> callback)
{
  pre_advance_callback_ = std::move(callback);
}

void
TimeDependentSourceSolver::SetPreAdvanceCallback(std::nullptr_t)
{
  pre_advance_callback_ = nullptr;
}

void
TimeDependentSourceSolver::SetPostAdvanceCallback(std::function<void()> callback)
{
  post_advance_callback_ = std::move(callback);
}

void
TimeDependentSourceSolver::SetPostAdvanceCallback(std::nullptr_t)
{
  post_advance_callback_ = nullptr;
}

} // namespace opensn
