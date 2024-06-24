// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_transient_solver/lbts_transient_solver.h"

#if 0
#include "modules/linear_boltzmann_solvers/lbs_solver/SourceFunctions/transient_source_function.h"
#include "modules/linear_boltzmann_solvers/B_DO_Solver/IterativeMethods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_linear_solver.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace lbs
{

DiscOrdTransientSolver::DiscOrdTransientSolver(const std::string& in_text_name)
  : DiscOrdKEigenvalueSolver(in_text_name)
{
  opensn::log.Log() << TextName() << " created.";
}

DiscOrdTransientSolver::~DiscOrdTransientSolver()
{
}

void
DiscOrdTransientSolver::Initialize()
{
  opensn::log.Log() << "Initializing " << TextName() << ".";
  options_.save_angular_flux = true;
  DiscOrdKEigenvalueSolver::Initialize();
  DiscOrdKEigenvalueSolver::Execute();

  if (transient_options_.verbosity_level >= 1)
  {
    const double FR = ComputeFissionRate(phi_new_local_);
    char buff[200];
    snprintf(buff, 200, " Initial Fission Rate FR=%12.6g", FR);
    opensn::log.Log() << TextName() << buff;
  }

  // Compute auxiliary vectors
  fission_rate_local_.resize(grid_ptr_->local_cells.size(), 0.0);
  phi_prev_local_ = phi_old_local_;
  precursor_prev_local_ = precursor_new_local_;
  psi_prev_local_ = psi_new_local_;

  if (transient_options_.verbosity_level >= 0)
  {
    const double beta = ComputeBeta();
    char buff[200];
    snprintf(buff,
             200,
             " Beta=%.2f [pcm] reactivity=%.3f [$]",
             beta * 1e5,
             (1.0 - 1.0 / GetKeff()) / beta);
    opensn::log.Log() << TextName() << buff;
  }

  // Initialize source func
  auto src_function = std::make_shared<TransientSourceFunction>(*this, this->dt_, this->method);

  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4, _5);
}

void
DiscOrdTransientSolver::Execute()
{
  opensn::log.Log() << "Executing " << TextName() << ".";

  const int max_num_steps = transient_options_.max_time_steps;
  const double max_time = transient_options_.t_final;
  int step_number = 0;
  while (((max_num_steps > 0 and step_number < max_num_steps) or (max_num_steps < 0)) and
         (time_ < max_time))
  {
    Step();

    PostStepCallBackFunction();

    if (not transient_options_.inhibit_advance)
    {
      Advance(); // new copied to prev + time+=dt
      ++step_number;
      transient_options_.inhibit_advance = false;
    }
  }

  UpdateFieldFunctions();

  opensn::log.Log() << "Done Executing " << TextName() << ".";
}

void
DiscOrdTransientSolver::Step()
{
  if (transient_options_.verbosity_level >= 2)
    opensn::log.Log() << TextName() << " Stepping with dt " << dt_;

  phi_old_local_ = phi_prev_local_;

  for (auto& groupset : groupsets_)
  {
    // Converge the scattering source with a fixed fission source and temporal source
    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    auto sweep_chunk = SetTransientSweepChunk(groupset);

    auto sweep_wgs_context_ptr = std::make_shared<SweepWGSContext<Mat, Vec, KSP>>(
      *this,
      groupset,
      active_set_source_function_,
      APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,
      APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES | APPLY_AGS_FISSION_SOURCES,
      options_.verbose_inner_iterations,
      sweep_chunk);

    WGSLinearSolver<Mat, Vec, KSP> solver(sweep_wgs_context_ptr);
    solver.Setup();
    solver.Solve();

     opensn::mpi_comm.barrier();
  }

  // Compute t^{n+1} value
  {
    const auto& BackwardEuler = opensn::SteppingMethod::IMPLICIT_EULER;
    const auto& CrankNicolson = opensn::SteppingMethod::CRANK_NICOLSON;

    double theta;
    if (method == BackwardEuler) theta = 1.0;
    else if (method == CrankNicolson)
      theta = 0.5;
    else
      theta = 0.7;
    const double inv_theta = 1.0 / theta;

    auto& phi = phi_new_local_;
    const auto& phi_prev = phi_prev_local_;
    for (size_t i = 0; i < phi.size(); ++i)
      phi[i] = inv_theta * (phi[i] + (theta - 1.0) * phi_prev[i]);

    if (options_.use_precursors) StepPrecursors();
  }

  const double FR_new = ComputeFissionProduction(phi_new_local_);

  // Print end of timestep
  if (transient_options_.verbosity_level >= 1)
  {
    char buff[200];
    snprintf(buff, 200, " dt=%.1e time=%10.4g FR=%12.6g", dt_, time_ + dt_, FR_new);
    opensn::log.Log() << TextName() << buff;
  }

  UpdateFieldFunctions();
}

void
DiscOrdTransientSolver::Advance()
{
  time_ += dt_;
  phi_prev_local_ = phi_new_local_;
  psi_prev_local_ = psi_new_local_;
  if (options_.use_precursors) precursor_prev_local_ = precursor_new_local_;
}

} // namespace lbs
#endif
