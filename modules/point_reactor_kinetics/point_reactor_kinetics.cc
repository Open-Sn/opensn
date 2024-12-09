// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/point_reactor_kinetics/point_reactor_kinetics.h"
#include "framework/physics/time_steppers/time_stepper.h"
#include "framework/event_system/physics_event_publisher.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/math/math.h"
#include <numeric>

namespace opensn
{

OpenSnRegisterObjectInNamespace(prk, PRKSolver);

InputParameters
PRKSolver::GetInputParameters()
{
  InputParameters params = opensn::Solver::GetInputParameters();

  params.SetDocGroup("prk");

  params.ChangeExistingParamToOptional("name", "PRKSolver");

  std::vector<double> default_lambdas = {0.0124, 0.0304, 0.111, 0.301, 1.14, 3.01};
  std::vector<double> default_betas = {0.00021, 0.00142, 0.00127, 0.00257, 0.00075, 0.00027};

  params.AddOptionalParameterArray(
    "precursor_lambdas", default_lambdas, "An array of decay constants");
  params.AddOptionalParameterArray(
    "precursor_betas", default_betas, "An array of fractional delayed neutron fractions");

  params.AddOptionalParameter("gen_time", 1.0e-5, "Neutron generation time [s]");
  params.AddOptionalParameter("initial_rho", 0.0, "Initial reactivity [$]");
  params.AddOptionalParameter("initial_source", 1.0, "Initial source strength [/s]");

  params.AddOptionalParameter("initial_population", 1.0, "Initial neutron population");

  params.AddOptionalParameter(
    "time_integration", "implicit_euler", "Time integration scheme to use");

  auto time_intgl_list =
    AllowableRangeList::New({"explicit_euler", "implicit_euler", "crank_nicolson"});

  params.ConstrainParameterRange("time_integration", std::move(time_intgl_list));

  params.ConstrainParameterRange("gen_time", AllowableRangeLowLimit::New(1.0e-12));
  params.ConstrainParameterRange("initial_source", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("initial_population", AllowableRangeLowLimit::New(0.0));
  return params;
}

std::shared_ptr<PRKSolver>
PRKSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  auto obj = factory.Create<PRKSolver>("prk::PRKSolver", params);
  return obj;
}

PRKSolver::PRKSolver(const InputParameters& params)
  : opensn::Solver(params),
    lambdas_(params.GetParamVectorValue<double>("precursor_lambdas")),
    betas_(params.GetParamVectorValue<double>("precursor_betas")),
    gen_time_(params.GetParamValue<double>("gen_time")),
    rho_(params.GetParamValue<double>("initial_rho")),
    source_strength_(params.GetParamValue<double>("initial_source")),
    time_integration_(params.GetParamValue<std::string>("time_integration")),
    num_precursors_(lambdas_.size())
{
  log.Log() << "Created solver " << Name();
  {
    std::stringstream outstr;
    outstr << "lambdas = ";
    for (double val : lambdas_)
      outstr << val << " ";
    log.Log() << outstr.str();
  }
  {
    std::stringstream outstr;
    outstr << "betas = ";
    for (double val : betas_)
      outstr << val << " ";
    log.Log() << outstr.str();
  }
}

void
PRKSolver::Initialize()
{
  // Check size
  if (lambdas_.size() != betas_.size())
    throw std::logic_error(Name() + ": Number of precursors cannot be "
                                    "deduced from precursor data because "
                                    "the data lists are of different size.");

  beta_ = std::accumulate(betas_.begin(), betas_.end(), 0.0);

  // Initializing linalg items
  const auto& J = num_precursors_;
  A_ = DenseMatrix<double>(J + 1, J + 1, 0.);
  I_ = A_;
  I_.SetDiagonal(1.0);

  x_t_ = Vector<double>(J + 1, 0.);

  // Assembling system
  A_(0, 0) = beta_ * (rho_ - 1.0) / gen_time_;
  for (size_t j = 1; j <= J; ++j)
  {
    A_(0, j) = lambdas_[j - 1];
    A_(j, j) = -lambdas_[j - 1];
    A_(j, 0) = betas_[j - 1] / gen_time_;
  }

  q_ = Vector<double>(J + 1, 0.);
  q_(0) = source_strength_;

  // Initializing x
  // If there is a source and the reactivity is < 0 then
  // there exists a unique solution.
  if (source_strength_ > 0.0 and rho_ < 0.0)
  {
    const auto b_theta = Scaled(q_, -1.0);

    x_t_ = Mult(Inverse(A_), b_theta);
  }
  // Otherwise we initialize the system as a critical system with
  // no source.
  else
  {
    auto A_temp = A_;
    auto b_temp = x_t_;
    b_temp.Set(0.0);
    b_temp(0) = 1.0;
    for (unsigned int i = 0; i < A_temp.Columns(); ++i)
      A_temp(0, i) = 0.0;
    A_temp(0, 0) = 1.0;
    x_t_ = Mult(Inverse(A_temp), b_temp);
  }

  log.Log() << "Final: " << x_t_.PrintStr();
}

void
PRKSolver::Execute()
{
  auto& physics_ev_pub = PhysicsEventPublisher::GetInstance();

  while (timestepper_->IsActive())
  {
    physics_ev_pub.SolverStep(*this);
    physics_ev_pub.SolverAdvance(*this);
  }
}

void
PRKSolver::Step()
{
  log.Log() << "Solver \"" + Name() + "\" " + timestepper_->StringTimeInfo();

  const double dt = timestepper_->TimeStepSize();

  A_(0, 0) = beta_ * (rho_ - 1.0) / gen_time_;

  if (time_integration_ == "implicit_euler" or time_integration_ == "crank_nicolson")
  {
    double theta = 1.0;

    if (time_integration_ == "implicit_euler")
      theta = 1.0;
    else if (time_integration_ == "crank_nicolson")
      theta = 0.5;

    const double inv_tau = theta * dt;

    auto A_theta = Subtract(I_, Scaled(A_, inv_tau));
    auto b_theta = Add(x_t_, Scaled(q_, inv_tau));

    auto x_theta = Mult(Inverse(A_theta), b_theta);

    x_tp1_ = Add(x_t_, Scaled(Subtract(x_theta, x_t_), 1.0 / theta));
  }
  else if (time_integration_ == "explicit_euler")
  {
    x_tp1_ = Add(x_t_, Add(Scaled(Mult(A_, x_t_), dt), Scaled(q_, dt)));
  }
  else
    OpenSnLogicalError("Unsupported time integration scheme.");

  if ((std::abs(x_t_(0)) > 1e-12) && std::abs((x_tp1_(0) / x_t_(0)) - 1.) > 1e-12)
    period_tph_ = dt / std::log(x_tp1_(0) / x_t_(0));
  else
    period_tph_ = 0.0;

  if (period_tph_ > 0.0 and period_tph_ > 1.0e6)
    period_tph_ = 1.0e6;
  if (period_tph_ < 0.0 and period_tph_ < -1.0e6)
    period_tph_ = -1.0e6;
}

void
PRKSolver::Advance()
{
  x_t_ = x_tp1_;
  timestepper_->Advance();
}

ParameterBlock
PRKSolver::GetInfo(const ParameterBlock& params) const
{
  const auto param_name = params.GetParamValue<std::string>("name");

  if (param_name == "neutron_population")
    return ParameterBlock("", x_t_(0));
  else if (param_name == "population_next")
    return ParameterBlock("", PopulationNew());
  else if (param_name == "period")
    return ParameterBlock("", period_tph_);
  else if (param_name == "rho")
    return ParameterBlock("", rho_);
  else if (param_name == "solution")
  {
    std::vector<double> sln = x_t_.ToStdVector();
    return ParameterBlock("", sln);
  }
  else if (param_name == "time_integration")
    return ParameterBlock("", time_integration_);
  else if (param_name == "time_next")
    return ParameterBlock("", TimeNew());
  else if (param_name == "test_arb_info")
  {
    ParameterBlock block;

    block.AddParameter("name", Name());
    block.AddParameter("time_integration", time_integration_);
    block.AddParameter("rho", rho_);
    block.AddParameter("max_timesteps", timestepper_->MaxTimeSteps());

    return block;
  }
  else
    OpenSnInvalidArgument("Unsupported info name \"" + param_name + "\".");
}

double
PRKSolver::PopulationPrev() const
{
  return x_t_(0);
}

double
PRKSolver::PopulationNew() const
{
  return x_tp1_(0);
}

double
PRKSolver::Period() const
{
  return period_tph_;
}

double
PRKSolver::TimePrev() const
{
  return timestepper_->Time();
}

double
PRKSolver::TimeNew() const
{
  return timestepper_->Time() + timestepper_->TimeStepSize();
}

Vector<double>
PRKSolver::SolutionPrev() const
{
  return x_t_;
}

Vector<double>
PRKSolver::SolutionNew() const
{
  return x_tp1_;
}

void
PRKSolver::SetRho(double value)
{
  rho_ = value;
}

void
PRKSolver::SetProperties(const ParameterBlock& params)
{
  opensn::Solver::SetProperties(params);

  for (const auto& param : params)
  {
    const std::string& param_name = param.GetName();
    if (param_name == "rho")
      SetRho(param.GetValue<double>());
  }
}

} // namespace opensn
