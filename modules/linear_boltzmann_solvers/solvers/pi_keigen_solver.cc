// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT
#include "framework/logging/log_exceptions.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/utils/hdf_utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/discrete_ordinates_keigen_acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/ags_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include <iomanip>
#include "sys/stat.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, PowerIterationKEigenSolver);

InputParameters
PowerIterationKEigenSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();

  params.SetGeneralDescription("Implementation of a k-Eigenvalue solver using Power Iteration");
  params.ChangeExistingParamToOptional("name", "PowerIterationKEigenSolver");
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "An existing discrete ordinates problem");
  params.AddOptionalParameter<std::shared_ptr<DiscreteOrdinatesKEigenAcceleration>>(
    "acceleration", {}, "The acceleration method");
  params.AddOptionalParameter("max_iters", 1000, "Maximum power iterations allowed");
  params.AddOptionalParameter("k_tol", 1.0e-10, "Tolerance on the k-eigenvalue");
  params.AddOptionalParameter(
    "reset_solution", true, "If set to true will initialize the flux moments to 1.0");
  params.AddOptionalParameter("reset_phi0", true, "If true, reinitializes scalar fluxes to 1.0");
  return params;
}

std::shared_ptr<PowerIterationKEigenSolver>
PowerIterationKEigenSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<PowerIterationKEigenSolver>("lbs::PowerIterationKEigenSolver", params);
}

PowerIterationKEigenSolver::PowerIterationKEigenSolver(const InputParameters& params)
  : Solver(params),
    do_problem_(params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>(("problem"))),
    acceleration_(
      params.GetSharedPtrParam<DiscreteOrdinatesKEigenAcceleration>("acceleration", false)),
    max_iters_(params.GetParamValue<size_t>("max_iters")),
    k_eff_(1.0),
    k_tolerance_(params.GetParamValue<double>("k_tol")),
    F_prev_(1.0),
    reset_phi0_(params.GetParamValue<bool>("reset_phi0")),
    q_moments_local_(do_problem_->GetQMomentsLocal()),
    phi_old_local_(do_problem_->GetPhiOldLocal()),
    phi_new_local_(do_problem_->GetPhiNewLocal()),
    groupsets_(do_problem_->GetGroupsets()),
    front_gs_(groupsets_.front())
{
}

void
PowerIterationKEigenSolver::Initialize()
{
  do_problem_->Initialize();

  auto& options = do_problem_->GetOptions();
  active_set_source_function_ = do_problem_->GetActiveSetSourceFunction();
  ags_solver_ = do_problem_->GetAGSSolver();

  for (auto& wgs_solver : do_problem_->GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<WGSContext>(context);

    OpenSnLogicalErrorIf(not wgs_context, ": Cast failed");

    wgs_context->lhs_src_scope.Unset(APPLY_WGS_FISSION_SOURCES); // lhs_scope
    wgs_context->rhs_src_scope.Unset(APPLY_AGS_FISSION_SOURCES); // rhs_scope
  }

  ags_solver_->SetVerbosity(options.verbose_ags_iterations);

  front_wgs_solver_ = do_problem_->GetWGSSolvers().at(front_gs_.id);
  front_wgs_context_ = std::dynamic_pointer_cast<WGSContext>(front_wgs_solver_->GetContext());

  OpenSnLogicalErrorIf(not front_wgs_context_, ": Casting failed");

  bool restart_successful = false;
  if (not options.read_restart_path.empty())
    restart_successful = ReadRestartData();

  if (reset_phi0_ and not restart_successful)
    LBSVecOps::SetPhiVectorScalarValues(*do_problem_, PhiSTLOption::PHI_OLD, 1.0);

  F_prev_ = ComputeFissionProduction(*do_problem_, phi_old_local_);

  if (acceleration_)
    acceleration_->Initialize(*this);
}

void
PowerIterationKEigenSolver::Execute()
{
  if (acceleration_)
    acceleration_->PreExecute();

  auto& options = do_problem_->GetOptions();
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  // Start power iterations
  int nit = 0;
  bool converged = false;
  while (nit < max_iters_)
  {
    // Set the fission source
    SetLBSFissionSource(phi_old_local_, false);
    Scale(q_moments_local_, 1.0 / k_eff_);

    if (acceleration_)
      acceleration_->PrePowerIteration();

    // This solves the inners for transport
    ags_solver_->Solve();

    // Get k-eigenvalue from the acceleration method
    if (acceleration_)
    {
      k_eff_ = acceleration_->PostPowerIteration();
    }
    // Recompute k-eigenvalue
    else
    {
      const auto F_new = ComputeFissionProduction(*do_problem_, phi_new_local_);
      k_eff_ = F_new / F_prev_ * k_eff_;
      F_prev_ = F_new;
    }

    const double reactivity = (k_eff_ - 1.0) / k_eff_;

    // Check convergence, bookkeeping
    k_eff_change = fabs(k_eff_ - k_eff_prev) / k_eff_;
    k_eff_prev = k_eff_;
    nit += 1;

    converged = k_eff_change < std::max(k_tolerance_, 1.0e-12);

    // Print iteration summary
    if (options.verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info << program_timer.GetTimeString() << " "
                  << "  Iteration " << std::setw(5) << nit << "  k_eff " << std::setw(11)
                  << std::setprecision(7) << k_eff_ << "  k_eff change " << std::setw(12)
                  << k_eff_change << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged)
        k_iter_info << " CONVERGED\n";

      log.Log() << k_iter_info.str();
    }

    if (options.restart_writes_enabled and do_problem_->TriggerRestartDump())
      WriteRestartData();

    if (converged)
      break;
  } // for k iterations

  // If restarts are enabled, always write a restart dump upon convergence or
  // when we reach the iteration limit
  if (options.restart_writes_enabled)
    WriteRestartData();

  // Print summary
  int total_num_sweeps = 0;
  for (auto& wgs_solver : do_problem_->GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<WGSContext>(context);
    total_num_sweeps += wgs_context->counter_applications_of_inv_op;
  }

  log.Log() << "\n";
  log.Log() << "        Final k-eigenvalue    :        " << std::setprecision(7) << k_eff_;
  log.Log() << "        Final change          :        " << std::setprecision(6) << k_eff_change
            << " (Total number of sweeps:" << total_num_sweeps << ")"
            << "\n\n";

  if (options.use_precursors)
  {
    ComputePrecursors(*do_problem_);
    Scale(do_problem_->GetPrecursorsNewLocal(), 1.0 / k_eff_);
  }

  do_problem_->UpdateFieldFunctions();

  log.Log() << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}

void
PowerIterationKEigenSolver::SetLBSFissionSource(const std::vector<double>& input,
                                                const bool additive)
{
  if (not additive)
    Set(q_moments_local_, 0.0);

  for (auto& groupset : groupsets_)
  {
    active_set_source_function_(
      groupset, q_moments_local_, input, APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  }
}

void
PowerIterationKEigenSolver::SetLBSScatterSource(const std::vector<double>& input,
                                                const bool additive,
                                                const bool suppress_wg_scat)
{
  if (not additive)
    Set(q_moments_local_, 0.0);

  SourceFlags source_flags = APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES;
  if (suppress_wg_scat)
    source_flags |= SUPPRESS_WG_SCATTER;

  for (auto& groupset : groupsets_)
  {
    active_set_source_function_(groupset, q_moments_local_, input, source_flags);
  }
}

bool
PowerIterationKEigenSolver::ReadRestartData()
{
  auto& fname = do_problem_->GetOptions().read_restart_path;
  auto& phi_old_local = do_problem_->GetPhiOldLocal();
  auto& groupsets = do_problem_->GetGroupsets();

  auto file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  bool success = (file >= 0);
  if (file >= 0)
  {
    // Read phi
    success &= H5ReadDataset1D<double>(file, "phi_old", phi_old_local);

    // Read psi
    int gs_id = 0;
    for (auto gs : groupsets)
    {
      if (gs.angle_agg)
      {
        std::string name = "delayed_psi_old_gs" + std::to_string(gs_id);
        if (H5Has(file, name))
        {
          std::vector<double> psi;
          success &= H5ReadDataset1D<double>(file, name.c_str(), psi);
          gs.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(psi);
        }
      }
      ++gs_id;
    }

    // Read keff and Fprev
    success &= H5ReadAttribute<double>(file, "keff", k_eff_) and
               H5ReadAttribute<double>(file, "Fprev", F_prev_);

    H5Fclose(file);
  }

  if (success)
    log.Log() << "Successfully read restart data." << std::endl;
  else
    log.Log() << "Failed to read restart data." << std::endl;

  return success;
}

bool
PowerIterationKEigenSolver::WriteRestartData()
{
  auto& options = do_problem_->GetOptions();
  auto fname = options.write_restart_path;
  auto& phi_old_local = do_problem_->GetPhiOldLocal();
  auto& groupsets = do_problem_->GetGroupsets();

  auto file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  bool success = (file >= 0);
  if (file >= 0)
  {
    // Write phi
    success &= H5WriteDataset1D<double>(file, "phi_old", phi_old_local);

    // Write psi
    if (options.write_delayed_psi_to_restart)
    {
      int gs_id = 0;
      for (auto gs : do_problem_->GetGroupsets())
      {
        if (gs.angle_agg)
        {
          auto psi = gs.angle_agg->GetOldDelayedAngularDOFsAsSTLVector();
          if (not psi.empty())
          {
            std::string name = "delayed_psi_old_gs" + std::to_string(gs_id);
            success &= H5WriteDataset1D<double>(file, name, psi);
          }
        }
        ++gs_id;
      }
    }

    success &= H5CreateAttribute<double>(file, "keff", k_eff_) and
               H5CreateAttribute<double>(file, "Fprev", F_prev_);

    H5Fclose(file);
  }

  if (success)
  {
    do_problem_->UpdateRestartWriteTime();
    log.Log() << "Successfully wrote restart data." << std::endl;
  }
  else
    log.Log() << "Failed to write restart data." << std::endl;

  return success;
}

} // namespace opensn
