// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/preconditioning/lbs_shell_operations.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <petscksp.h>
#include "caliper/cali.h"
#include <iomanip>

namespace opensn
{

using PCShellPtr = PetscErrorCode (*)(PC, Vec, Vec);
using namespace std::chrono;

SweepWGSContext::SweepWGSContext(DiscreteOrdinatesProblem& do_problem,
                                 LBSGroupset& groupset,
                                 const SetSourceFunction& set_source_function,
                                 SourceFlags lhs_scope,
                                 SourceFlags rhs_scope,
                                 bool log_info,
                                 std::shared_ptr<SweepChunk> swp_chnk)
  : WGSContext(do_problem, groupset, set_source_function, lhs_scope, rhs_scope, log_info),
    sweep_chunk(std::move(swp_chnk)),
    sweep_scheduler(do_problem.GetSweepType() == "AAH" ? SchedulingAlgorithm::DEPTH_OF_GRAPH
                                                       : SchedulingAlgorithm::FIRST_IN_FIRST_OUT,
                    *groupset.angle_agg,
                    *sweep_chunk)
{
}

void
SweepWGSContext::PreSetupCallback()
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::PreSetupCallback");

  if (log_info)
    log.Log() << "\n\n"
              << "********** Solving groupset " << groupset.id << " with "
              << LinearSolver::IterativeMethodName(groupset.iterative_method) << ".\n\n"
              << "Quadrature number of angles: " << groupset.quadrature->abscissae.size() << "\n"
              << "Groups " << groupset.groups.front().id << " " << groupset.groups.back().id
              << "\n\n";
}

void
SweepWGSContext::SetPreconditioner(KSP& solver)
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::SetPreconditioner");

  auto& ksp = solver;

  PC pc;
  KSPGetPC(ksp, &pc);

  if (groupset.apply_wgdsa or groupset.apply_tgdsa)
  {
    PCSetType(pc, PCSHELL);
    PCShellSetApply(pc, (PCShellPtr)WGDSA_TGDSA_PreConditionerMult);
    PCShellSetContext(pc, &(*this));
  }

  KSPSetPCSide(ksp, PC_LEFT);
  KSPSetUp(ksp);
}

std::pair<int64_t, int64_t>
SweepWGSContext::GetSystemSize()
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::SystemSize");

  const size_t local_node_count = do_problem.GetLocalNodeCount();
  const size_t global_node_count = do_problem.GetGlobalNodeCount();
  const size_t num_moments = do_problem.GetNumMoments();

  const size_t groupset_numgrps = groupset.groups.size();
  const auto num_delayed_psi_info = groupset.angle_agg->GetNumDelayedAngularDOFs();
  const size_t local_size =
    local_node_count * num_moments * groupset_numgrps + num_delayed_psi_info.first;
  const size_t global_size =
    global_node_count * num_moments * groupset_numgrps + num_delayed_psi_info.second;
  const size_t num_angles = groupset.quadrature->abscissae.size();
  const size_t num_psi_global = global_node_count * num_angles * groupset.groups.size();
  const size_t num_delayed_psi_global = num_delayed_psi_info.second;

  if (log_info)
  {
    log.Log() << "Total number of angular unknowns: " << num_psi_global << "\n"
              << "Number of lagged angular unknowns: " << num_delayed_psi_global << "("
              << std::setprecision(2)
              << static_cast<double>(num_delayed_psi_global) * 100 /
                   static_cast<double>(num_psi_global)
              << "%)";
  }

  return {static_cast<int64_t>(local_size), static_cast<int64_t>(global_size)};
}

void
SweepWGSContext::ApplyInverseTransportOperator(SourceFlags scope)
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::ApplyInverseTransportOperator");

  ++counter_applications_of_inv_op;

  // Sweep
  const bool use_bndry_source_flag =
    (scope & APPLY_FIXED_SOURCES) and (not do_problem.GetOptions().use_src_moments);
  const bool zero_incoming_delayed_psi = (scope & ZERO_INCOMING_DELAYED_PSI);
  dynamic_cast<DiscreteOrdinatesProblem&>(do_problem).ZeroOutflowBalanceVars(groupset);
  sweep_scheduler.PrepareForSweep(use_bndry_source_flag, zero_incoming_delayed_psi);

  high_resolution_clock::time_point sweep_start = high_resolution_clock::now();
  sweep_scheduler.Sweep();
  high_resolution_clock::time_point sweep_end = high_resolution_clock::now();

  double sweep_time = (duration_cast<nanoseconds>(sweep_end - sweep_start).count()) / 1.0e+9;
  sweep_times.push_back(sweep_time);
}

void
SweepWGSContext::PostSolveCallback()
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::PostSolveCallback");

  // Perform final sweep with converged phi and delayed psi dofs. This step is necessary for
  // Krylov methods to recover the actual solution (this includes all of the PETSc methods
  // currently used in OpenSn).
  if (groupset.iterative_method == LinearSolver::IterativeMethod::PETSC_GMRES or
      groupset.iterative_method == LinearSolver::IterativeMethod::PETSC_BICGSTAB or
      (groupset.iterative_method == LinearSolver::IterativeMethod::PETSC_RICHARDSON and
       groupset.max_iterations > 1))
  {
    const auto scope = lhs_src_scope | rhs_src_scope;
    set_source_function(
      groupset, do_problem.GetQMomentsLocal(), do_problem.GetPhiOldLocal(), scope);
    ApplyInverseTransportOperator(scope);
    LBSVecOps::GSScopedCopyPrimarySTLvectors(
      do_problem, groupset, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);
  }

  if (log_info)
  {
    double tot_sweep_time = 0.0;
    auto num_sweeps = static_cast<double>(sweep_times.size());
    for (auto time : sweep_times)
      tot_sweep_time += time;
    double avg_sweep_time = tot_sweep_time / num_sweeps;
    size_t num_angles = groupset.quadrature->abscissae.size();
    size_t num_unknowns = do_problem.GetGlobalNodeCount() * num_angles * groupset.groups.size();

    log.Log() << "\n       Average sweep time (s):        "
              << tot_sweep_time / static_cast<double>(sweep_times.size())
              << "\n       Sweep Time/Unknown (ns):       "
              << avg_sweep_time * 1.0e9 * opensn::mpi_comm.size() /
                   static_cast<double>(num_unknowns)
              << "\n       Number of unknowns per sweep:  " << num_unknowns << "\n\n";
  }
}

} // namespace opensn
