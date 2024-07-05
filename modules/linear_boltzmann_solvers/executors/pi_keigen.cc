// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/executors/pi_keigen.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_linear_solver.h"
#include "framework/logging/log_exceptions.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/utils/hdf_utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include <iomanip>
#include "sys/stat.h"

namespace opensn
{
namespace lbs
{

OpenSnRegisterObjectInNamespace(lbs, PowerIterationKEigen);

InputParameters
PowerIterationKEigen::GetInputParameters()
{
  InputParameters params = opensn::Solver::GetInputParameters();

  params.SetGeneralDescription("Implementation of a k-Eigenvalue solver using Power Iteration");
  params.SetDocGroup("LBSExecutors");
  params.ChangeExistingParamToOptional("name", "PowerIterationKEigen");
  params.AddRequiredParameter<size_t>("lbs_solver_handle", "Handle to an existing lbs solver");
  params.AddOptionalParameter("max_iters", 1000, "Maximum power iterations allowed");
  params.AddOptionalParameter("k_tol", 1.0e-10, "Tolerance on the k-eigenvalue");
  params.AddOptionalParameter(
    "reset_solution", true, "If set to true will initialize the flux moments to 1.0");
  params.AddOptionalParameter("reset_phi0", true, "If true, reinitializes scalar fluxes to 1.0");

  return params;
}

PowerIterationKEigen::PowerIterationKEigen(const InputParameters& params)
  : opensn::Solver(params),
    lbs_solver_(
      GetStackItem<LBSSolver>(object_stack, params.GetParamValue<size_t>("lbs_solver_handle"))),
    max_iters_(params.GetParamValue<size_t>("max_iters")),
    k_eff_(1.0),
    k_tolerance_(params.GetParamValue<double>("k_tol")),
    F_prev_(1.0),
    reset_phi0_(params.GetParamValue<bool>("reset_phi0")),
    q_moments_local_(lbs_solver_.QMomentsLocal()),
    phi_old_local_(lbs_solver_.PhiOldLocal()),
    phi_new_local_(lbs_solver_.PhiNewLocal()),
    groupsets_(lbs_solver_.Groupsets()),
    front_gs_(groupsets_.front())
{
  lbs_solver_.Options().enable_ags_restart_write = false;
}

void
PowerIterationKEigen::Initialize()
{
  lbs_solver_.Initialize();

  active_set_source_function_ = lbs_solver_.GetActiveSetSourceFunction();
  primary_ags_solver_ = lbs_solver_.GetPrimaryAGSSolver();

  for (auto& wgs_solver : lbs_solver_.GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<lbs::WGSContext>(context);

    OpenSnLogicalErrorIf(not wgs_context, ": Cast failed");

    wgs_context->lhs_src_scope_.Unset(APPLY_WGS_FISSION_SOURCES); // lhs_scope
    wgs_context->rhs_src_scope_.Unset(APPLY_AGS_FISSION_SOURCES); // rhs_scope
  }

  primary_ags_solver_->SetVerbosity(lbs_solver_.Options().verbose_ags_iterations);

  front_wgs_solver_ = lbs_solver_.GetWGSSolvers().at(front_gs_.id_);
  front_wgs_context_ = std::dynamic_pointer_cast<lbs::WGSContext>(front_wgs_solver_->GetContext());

  OpenSnLogicalErrorIf(not front_wgs_context_, ": Casting failed");

  if (reset_phi0_ and lbs_solver_.Options().read_restart_path.empty())
    lbs_solver_.SetPhiVectorScalarValues(phi_old_local_, 1.0);
}

void
PowerIterationKEigen::Execute()
{
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  if (not lbs_solver_.Options().read_restart_path.empty())
    ReadRestartData();

  // Start power iterations
  int nit = 0;
  bool converged = false;
  while (nit < max_iters_)
  {
    // Set the fission source
    SetLBSFissionSource(phi_old_local_, false);
    Scale(q_moments_local_, 1.0 / k_eff_);

    // This solves the inners for transport
    primary_ags_solver_->Setup();
    primary_ags_solver_->Solve();

    // Recompute k-eigenvalue
    double F_new = lbs_solver_.ComputeFissionProduction(phi_new_local_);
    k_eff_ = F_new / F_prev_ * k_eff_;
    double reactivity = (k_eff_ - 1.0) / k_eff_;

    // Check convergence, bookkeeping
    k_eff_change = fabs(k_eff_ - k_eff_prev) / k_eff_;
    k_eff_prev = k_eff_;
    F_prev_ = F_new;
    nit += 1;

    if (k_eff_change < std::max(k_tolerance_, 1.0e-12))
      converged = true;

    // Print iteration summary
    if (lbs_solver_.Options().verbose_outer_iterations)
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

    if (lbs_solver_.RestartsEnabled() and lbs_solver_.TriggerRestartDump())
      WriteRestartData();

    if (converged)
      break;
  } // for k iterations

  // If restarts are enabled, always write a restart dump upon convergence or
  // when we reach the iteration limit
  if (lbs_solver_.RestartsEnabled())
    WriteRestartData();

  // Print summary
  int total_num_sweeps = 0;
  for (auto& wgs_solver : lbs_solver_.GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<lbs::WGSContext>(context);
    total_num_sweeps += wgs_context->counter_applications_of_inv_op_;
  }

  log.Log() << "\n";
  log.Log() << "        Final k-eigenvalue    :        " << std::setprecision(7) << k_eff_;
  log.Log() << "        Final change          :        " << std::setprecision(6) << k_eff_change
            << " (Total number of sweeps:" << total_num_sweeps << ")"
            << "\n\n";

  if (lbs_solver_.Options().use_precursors)
  {
    lbs_solver_.ComputePrecursors();
    Scale(lbs_solver_.PrecursorsNewLocal(), 1.0 / k_eff_);
  }

  lbs_solver_.UpdateFieldFunctions();

  log.Log() << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}

void
PowerIterationKEigen::SetLBSFissionSource(const std::vector<double>& input, const bool additive)
{
  if (not additive)
    Set(q_moments_local_, 0.0);

  for (auto& groupset : groupsets_)
  {
    active_set_source_function_(groupset,
                                q_moments_local_,
                                input,
                                lbs_solver_.DensitiesLocal(),
                                APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  }
}

void
PowerIterationKEigen::SetLBSScatterSource(const std::vector<double>& input,
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
    active_set_source_function_(
      groupset, q_moments_local_, input, lbs_solver_.DensitiesLocal(), source_flags);
  }
}

void
PowerIterationKEigen::WriteRestartData()
{
  std::string fbase = lbs_solver_.Options().write_restart_path.string();
  std::string fname = fbase + std::to_string(opensn::mpi_comm.rank()) + ".restart.h5";

  // Write data
  bool location_succeeded = true;
  hid_t file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (not file)
    location_succeeded = false;
  else
  {
    location_succeeded = H5WriteDataset1D<double>(file, "phi_old", lbs_solver_.PhiOldLocal()) and
                         H5CreateAttribute<double>(file, "keff", k_eff_) and
                         H5CreateAttribute<double>(file, "Fprev", F_prev_);
    H5Fclose(file);
  }

  bool global_succeeded = true;
  mpi_comm.all_reduce(location_succeeded, global_succeeded, mpi::op::logical_and<bool>());
  if (global_succeeded)
  {
    log.Log() << "Successfully wrote restart data to " << fbase << "X.restart.h5";
    lbs_solver_.UpdateLastRestartWriteTime();
  }
  else
    log.Log0Error() << "Failed to write restart data to " << fbase << "X.restart.h5";
}

void
PowerIterationKEigen::ReadRestartData()
{
  std::string fbase = lbs_solver_.Options().read_restart_path.string();
  std::string fname = fbase + std::to_string(opensn::mpi_comm.rank()) + ".restart.h5";

  bool location_succeeded = true;
  hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (not file)
    location_succeeded = false;
  else
  {
    auto& phi_old_local = lbs_solver_.PhiOldLocal();
    phi_old_local.clear();
    phi_old_local = H5ReadDataset1D<double>(file, "phi_old");
    location_succeeded = (not phi_old_local.empty()) and
                         H5ReadAttribute<double>(file, "keff", k_eff_) and
                         H5ReadAttribute<double>(file, "Fprev", F_prev_);
    H5Fclose(file);
  }

  bool global_succeeded = true;
  mpi_comm.all_reduce(location_succeeded, global_succeeded, mpi::op::logical_and<bool>());
  if (global_succeeded)
    log.Log() << "Successfully read restart data from " << fbase << "X.restart.h5";
  else
    throw std::invalid_argument("Failed to read restart data from " + fbase + "X.restart.h5");
}

} // namespace lbs
} // namespace opensn
