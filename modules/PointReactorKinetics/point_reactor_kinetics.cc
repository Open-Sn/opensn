#include "opensn/modules/PointReactorKinetics/point_reactor_kinetics.h"

#include "opensn/framework/physics/TimeSteppers/TimeStepper.h"
#include "opensn/framework/physics/PhysicsEventPublisher.h"

#include "opensn/framework/ChiObjectFactory.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"

#include <numeric>

namespace prk
{

RegisterChiObject(prk, TransientSolver);

/**Sets input parameters.*/
chi::InputParameters
TransientSolver::GetInputParameters()
{
  chi::InputParameters params = chi_physics::Solver::GetInputParameters();

  params.SetDocGroup("prk");

  params.ChangeExistingParamToOptional("name", "prk_TransientSolver");

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

  using namespace chi_data_types;
  auto time_intgl_list =
    AllowableRangeList::New({"explicit_euler", "implicit_euler", "crank_nicolson"});

  params.ConstrainParameterRange("time_integration", std::move(time_intgl_list));

  params.ConstrainParameterRange("gen_time", AllowableRangeLowLimit::New(1.0e-12));
  params.ConstrainParameterRange("initial_source", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("initial_population", AllowableRangeLowLimit::New(0.0));
  return params;
}

/**Constructor.*/
TransientSolver::TransientSolver(const chi::InputParameters& params)
  : chi_physics::Solver(params.GetParamValue<std::string>("name")),
    lambdas_(params.GetParamVectorValue<double>("precursor_lambdas")),
    betas_(params.GetParamVectorValue<double>("precursor_betas")),
    gen_time_(params.GetParamValue<double>("gen_time")),
    rho_(params.GetParamValue<double>("initial_rho")),
    source_strength_(params.GetParamValue<double>("initial_source")),
    time_integration_(params.GetParamValue<std::string>("time_integration")),
    num_precursors_(lambdas_.size())
{
  Chi::log.Log() << "Created solver " << TextName();
  {
    std::stringstream outstr;
    outstr << "lambdas = ";
    for (double val : lambdas_)
      outstr << val << " ";
    Chi::log.Log() << outstr.str();
  }
  {
    std::stringstream outstr;
    outstr << "betas = ";
    for (double val : betas_)
      outstr << val << " ";
    Chi::log.Log() << outstr.str();
  }
}

/**Initialize function.*/
void
TransientSolver::Initialize()
{
  // Check size
  if (lambdas_.size() != betas_.size())
    throw std::logic_error(TextName() + ": Number of precursors cannot be "
                                        "deduced from precursor data because "
                                        "the data lists are of different size.");

  beta_ = std::accumulate(betas_.begin(), betas_.end(), /*init_val=*/0.0);

  // ============================= Initializing linalg items
  const auto& J = num_precursors_;
  A_ = chi_math::DynamicMatrix<double>(J + 1, J + 1, 0.0);
  I_ = A_;
  I_.SetDiagonal(1.0);

  x_t_ = chi_math::DynamicVector<double>(J + 1, 0.0);

  // ============================= Assembling system
  A_[0][0] = beta_ * (rho_ - 1.0) / gen_time_;
  for (size_t j = 1; j <= J; ++j)
  {
    A_[0][j] = lambdas_[j - 1];
    A_[j][j] = -lambdas_[j - 1];
    A_[j][0] = betas_[j - 1] / gen_time_;
  }

  q_.resize(J + 1, 0.0);
  q_[0] = source_strength_;

  // ============================= Initializing x
  // If there is a source and the reactivity is < 0 then
  // there exists a unique solution.
  if (source_strength_ > 0.0 and rho_ < 0.0)
  {
    const auto b_theta = -1.0 * q_;

    x_t_ = A_.Inverse() * b_theta;
  }
  // Otherwise we initialize the system as a critical system with
  // no source.
  else
  {
    auto A_temp = A_;
    auto b_temp = x_t_;
    b_temp.Set(0.0);
    b_temp[0] = 1.0;
    for (auto& val : A_temp[0])
      val = 0.0;
    A_temp[0][0] = 1.0;

    x_t_ = A_temp.Inverse() * b_temp;
  }

  Chi::log.Log() << "Final: " << x_t_.PrintStr();
}

/**Execution function.*/
void
TransientSolver::Execute()
{
  auto& physics_ev_pub = chi_physics::PhysicsEventPublisher::GetInstance();

  while (timestepper_->IsActive())
  {
    physics_ev_pub.SolverStep(*this);
    physics_ev_pub.SolverAdvance(*this);
  }
}

/**Step function*/
void
TransientSolver::Step()
{
  Chi::log.Log() << "Solver \"" + TextName() + "\" " + timestepper_->StringTimeInfo();

  const double dt = timestepper_->TimeStepSize();

  A_[0][0] = beta_ * (rho_ - 1.0) / gen_time_;

  if (time_integration_ == "implicit_euler" or time_integration_ == "crank_nicolson")
  {
    double theta = 1.0;

    if (time_integration_ == "implicit_euler") theta = 1.0;
    else if (time_integration_ == "crank_nicolson")
      theta = 0.5;

    const double inv_tau = theta * dt;

    auto A_theta = I_ - inv_tau * A_;
    auto b_theta = x_t_ + inv_tau * q_;

    auto x_theta = A_theta.Inverse() * b_theta;

    x_tp1_ = x_t_ + (1.0 / theta) * (x_theta - x_t_);
  }
  else if (time_integration_ == "explicit_euler") { x_tp1_ = x_t_ + dt * A_ * x_t_ + dt * q_; }
  else
    ChiLogicalError("Unsupported time integration scheme.");

  period_tph_ = dt / log(x_tp1_[0] / x_t_[0]);

  if (period_tph_ > 0.0 and period_tph_ > 1.0e6) period_tph_ = 1.0e6;
  if (period_tph_ < 0.0 and period_tph_ < -1.0e6) period_tph_ = -1.0e6;
}

/**Advance time values function.*/
void
TransientSolver::Advance()
{
  x_t_ = x_tp1_;
  timestepper_->Advance();
}

chi::ParameterBlock
TransientSolver::GetInfo(const chi::ParameterBlock& params) const
{
  const auto param_name = params.GetParamValue<std::string>("name");

  if (param_name == "neutron_population") return chi::ParameterBlock("", x_t_[0]);
  else if (param_name == "population_next")
    return chi::ParameterBlock("", PopulationNew());
  else if (param_name == "period")
    return chi::ParameterBlock("", period_tph_);
  else if (param_name == "rho")
    return chi::ParameterBlock("", rho_);
  else if (param_name == "solution")
    return chi::ParameterBlock("", x_t_.elements_);
  else if (param_name == "time_integration")
    return chi::ParameterBlock("", time_integration_);
  else if (param_name == "time_next")
    return chi::ParameterBlock("", TimeNew());
  else if (param_name == "test_arb_info")
  {
    chi::ParameterBlock block;

    block.AddParameter("name", TextName());
    block.AddParameter("time_integration", time_integration_);
    block.AddParameter("rho", rho_);
    block.AddParameter("max_timesteps", timestepper_->MaxTimeSteps());

    return block;
  }
  else
    ChiInvalidArgument("Unsupported info name \"" + param_name + "\".");
}

// ##################################################################
/**Returns the population at the previous time step.*/
double
TransientSolver::PopulationPrev() const
{
  return x_t_[0];
}
/**Returns the population at the next time step.*/
double
TransientSolver::PopulationNew() const
{
  return x_tp1_[0];
}

/**Returns the period computed for the last time step.*/
double
TransientSolver::Period() const
{
  return period_tph_;
}

/**Returns the time computed for the last time step.*/
double
TransientSolver::TimePrev() const
{
  return timestepper_->Time();
}

/**Returns the time computed for the next time step.*/
double
TransientSolver::TimeNew() const
{
  return timestepper_->Time() + timestepper_->TimeStepSize();
}

/**Returns the solution at the previous time step.*/
std::vector<double>
TransientSolver::SolutionPrev() const
{
  return x_t_.elements_;
}
/**Returns the solution at the next time step.*/
std::vector<double>
TransientSolver::SolutionNew() const
{
  return x_tp1_.elements_;
}

// ##################################################################
/**Sets the value of rho.*/
void
TransientSolver::SetRho(double value)
{
  rho_ = value;
}

/**\addtogroup prk
 *
 * \section Properties Properties that can be set
 * The following properties can be set via the lua call
 * `chi_lua::chiSolverSetProperties`
 * \copydoc prk::TransientSolver::SetProperties
 * */

/** PRK Transient solver settable properties:
* - `rho`, The current reactivity

Parents:
\copydoc chi_physics::Solver::SetProperties*/
void
TransientSolver::SetProperties(const chi::ParameterBlock& params)
{
  chi_physics::Solver::SetProperties(params);

  for (const auto& param : params)
  {
    const std::string& param_name = param.Name();
    if (param_name == "rho") SetRho(param.GetValue<double>());
  }
}

} // namespace prk
