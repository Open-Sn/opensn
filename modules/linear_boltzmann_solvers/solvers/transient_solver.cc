// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/transient_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/source_functions/transient_source_function.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "framework/logging/log.h"

namespace opensn
{
namespace
{
double
ThetaFromMethod(SteppingMethod method)
{
  if (method == SteppingMethod::IMPLICIT_EULER)
    return 1.0;
  if (method == SteppingMethod::CRANK_NICOLSON)
    return 0.5;
  return 0.7;
}
} // namespace

OpenSnRegisterObjectInNamespace(lbs, TransientKEigenSolver);

InputParameters
TransientKEigenSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();
  params.SetGeneralDescription("Generalized implementation of a transient k-eigenvalue solver.");
  params.ChangeExistingParamToOptional("name", "TransientKEigenSolver");
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem", "An existing lbs problem");
  params.AddOptionalParameter<std::shared_ptr<DiscreteOrdinatesKEigenAcceleration>>(
    "acceleration", {}, "The acceleration method");
  params.AddOptionalParameter("max_iters", 1000, "Maximum power iterations allowed");
  params.AddOptionalParameter("k_tol", 1.0e-10, "Tolerance on the k-eigenvalue");
  params.AddOptionalParameter(
    "reset_solution", true, "If set to true will initialize the flux moments to 1.0");
  params.AddOptionalParameter("reset_phi0", true, "If true, reinitializes scalar fluxes to 1.0");

  return params;
}

void
TransientKEigenSolver::SetTimeStep(double dt)
{
  if (dt < 0.0)
    throw std::runtime_error(GetName() + " dt must be non-negative");
  dt_ = dt;
}

void
TransientKEigenSolver::SetTheta(double theta)
{
  if (theta < 0.0 or theta > 1.0)
    throw std::runtime_error(GetName() + " theta must be between 0.0 and 1.0.");
  theta_ = theta;
}

std::shared_ptr<TransientKEigenSolver>
TransientKEigenSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<TransientKEigenSolver>("lbs::TransientKEigenSolver", params);
}

TransientKEigenSolver::TransientKEigenSolver(const InputParameters& params)
  : PowerIterationKEigenSolver(params),
    do_problem_(params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>("problem")),
    phi_new_local_(do_problem_->GetPhiNewLocal()),
    precursor_new_local_(do_problem_->GetPrecursorsNewLocal()),
    psi_new_local_(do_problem_->GetPsiNewLocal())
{
}

void
TransientKEigenSolver::Initialize()
{
  log.Log() << "Initializing " << GetName() << ".";
  do_problem_->GetOptions().save_angular_flux = true;
  do_problem_->SetTime(time_);
  do_problem_->SetTimeStep(dt_);
  do_problem_->SetTheta(ThetaFromMethod(method));
  PowerIterationKEigenSolver::Initialize();
  PowerIterationKEigenSolver::Execute();

  if (transient_options_.verbosity_level >= 1)
  {
    const double FR = ComputeFissionRate(*do_problem_, phi_new_local_);
    char buff[200];
    snprintf(buff,200, " Initial Fission Rate FR=%12.6g", FR);
    log.Log() << GetName() << buff;
  }

  // Compute auxiliary vectors
  fission_rate_local_.resize(do_problem_->GetGrid()->local_cells.size(), 0.0);
  phi_prev_local_ = phi_old_local_;
  precursor_prev_local_ = precursor_new_local_;
  psi_prev_local_ = psi_new_local_;
  do_problem_->UpdatePsiOld();

  if (transient_options_.verbosity_level >= 0)
  {
    /*
    const double beta = ComputeBeta();
    char buff[200];
    snprintf(buff,200, " Beta=%.2f [pcm] reactivity=%.3f [$]",
            beta*1e5, (1.0- 1.0 / GetKeff()) / beta);
    log.Log() << GetName() << buff;
    */
  }

  // Initialize source function
  auto src_function = std::make_shared<TransientSourceFunction>(*do_problem_, dt_, method);

  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4);
}

void
TransientKEigenSolver::Execute()
{
  log.Log() << "Executing " << GetName() << ".";

  const int max_num_steps = transient_options_.max_time_steps;
  const double max_time = transient_options_.t_final;
  int step_number = 0;
  while (((max_num_steps > 0 and step_number < max_num_steps) or
         (max_num_steps < 0)) and (time_ < max_time))
  {
    Step();

    //PostStepCallBackFunction();

    if (not transient_options_.inhibit_advance)
    {
      Advance(); // new copied to prev + time+=dt
      ++step_number;
      transient_options_.inhibit_advance = false;
    }
  }

  do_problem_->UpdateFieldFunctions();

  log.Log() << "Done Executing " << GetName() << ".";
}

void
TransientKEigenSolver::Step()
{
  auto& options = do_problem_->GetOptions();

  if (transient_options_.verbosity_level >= 2)
    log.Log() << GetName() << " Stepping with dt " << dt_;

  const double theta = ThetaFromMethod(method);
  do_problem_->SetTime(time_);
  do_problem_->SetTimeStep(dt_);
  do_problem_->SetTheta(theta);

  ags_solver_ = do_problem_->GetAGSSolver();
  ags_solver_->Solve();


  phi_old_local_ = phi_prev_local_;

  /*
  for (auto& groupset : groupsets_)
  {
    // Converge the scattering source with a fixed fission source and temporal source
    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    std::shared_ptr<SweepChunk> sweep_chunk = SetTransientSweepChunk(groupset);
    auto sweep_wgs_context_ptr = std::make_shared<SweepWGSContext>(
      *do_problem_,
      groupset,
      active_set_source_function_,
      APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,
      APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES | APPLY_AGS_FISSION_SOURCES,
      options.verbose_inner_iterations,
      sweep_chunk);
    WGSLinearSolver solver(sweep_wgs_context_ptr);
    
    solver.Setup();
    solver.Solve();

    mpi_comm.barrier();
  }
  */

  // Compute t^{n+1} value
  {
    const double inv_theta = 1.0/theta;

    auto& phi = phi_new_local_;
    const auto& phi_prev = phi_prev_local_;
    for (size_t i = 0; i < phi.size(); ++i)
      phi[i] = inv_theta*(phi[i] + (theta-1.0) * phi_prev[i]);

    if (options.use_precursors)
      StepPrecursors();
  }

  const double FR_new = ComputeFissionProduction(*do_problem_, phi_new_local_);

  // Print end of timestep
  if (transient_options_.verbosity_level >= 1)
  {
    char buff[200];
    snprintf(buff, 200, " dt=%.1e time=%10.4f FR=%12.6g", dt_, time_ + dt_, FR_new);
    log.Log() << GetName() << buff;
  }

  do_problem_->UpdateFieldFunctions();
}

void
TransientKEigenSolver::Advance()
{
  auto& options = do_problem_->GetOptions();
  time_ += dt_;
  do_problem_->SetTime(time_);
  phi_prev_local_ = phi_new_local_;
  psi_prev_local_ = psi_new_local_;
  do_problem_->UpdatePsiOld();
  if (options.use_precursors)
    precursor_prev_local_ = precursor_new_local_;
}

void
TransientKEigenSolver::StepPrecursors()
{
  const double theta = ThetaFromMethod(method);

  const double eff_dt = theta * dt_;

  precursor_new_local_.assign(precursor_new_local_.size(), 0.0);

  // Uses phi_new and precursor_prev_local to compute precursor_new_local(theta-flavor)
  auto& transport_views = do_problem_->GetCellTransportViews();
  for (const auto& cell : do_problem_->GetGrid()->local_cells)
  {
    const auto& fe_values = do_problem_->GetUnitCellMatrices()[cell.local_id];
    const double cell_volume = transport_views[cell.local_id].GetVolume();
    const auto& xs = do_problem_->GetBlockID2XSMap().at(cell.block_id);
    const auto& precursors = xs->GetPrecursors();
    const auto& nu_delayed_sigma_f = xs->GetNuDelayedSigmaF();

    // Compute delayed fission rate
    double delayed_fission = 0.0;
    for (int i = 0; i < transport_views[cell.local_id].GetNumNodes(); ++i)
    {
      const size_t uk_map = transport_views[cell.local_id].MapDOF(i, 0, 0);
      const double node_V_fraction = fe_values.intV_shapeI(i)/cell_volume;

      for (int g = 0; g < do_problem_->GetNumGroups(); ++g)
        delayed_fission += nu_delayed_sigma_f[g] *
                           phi_new_local_[uk_map + g] *
                           node_V_fraction;
    }

    // Loop over precursors
    const auto& max_precursors = do_problem_->GetMaxPrecursorsPerMaterial();
    for (unsigned int j = 0; j < xs->GetNumPrecursors(); ++j)
    {
      const size_t dof_map = cell.local_id * max_precursors + j;
      const auto& precursor = precursors[j];
      const double coeff = 1.0 / (1.0 + eff_dt * precursor.decay_constant);

      //contribute last time step precursors
      precursor_new_local_[dof_map] = coeff * precursor_prev_local_[dof_map];

      //contribute delayed fission production
      precursor_new_local_[dof_map] +=
        coeff * eff_dt * precursor.fractional_yield * delayed_fission;
    }
  }//for cell

  // Compute t^{n+1} value
  {
    auto& Cj = precursor_new_local_;
    const auto& Cj_prev = precursor_prev_local_;

    const double inv_theta = 1.0/theta;
    for (size_t i = 0; i < Cj.size(); ++i)
      Cj[i] = inv_theta * (Cj[i] + (theta - 1.0) * Cj_prev[i]);
  }
}

} // namespace opensn
