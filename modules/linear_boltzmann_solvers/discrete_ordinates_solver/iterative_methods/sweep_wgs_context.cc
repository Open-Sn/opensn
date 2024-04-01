#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/preconditioning/lbs_shell_operations.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <petscksp.h>
#include "caliper/cali.h"
#include <iomanip>

namespace opensn
{
namespace lbs
{

typedef PetscErrorCode (*PCShellPtr)(PC, Vec, Vec);

SweepWGSContext::SweepWGSContext(DiscreteOrdinatesSolver& lbs_solver,
                                 LBSGroupset& groupset,
                                 const SetSourceFunction& set_source_function,
                                 SourceFlags lhs_scope,
                                 SourceFlags rhs_scope,
                                 bool log_info,
                                 std::shared_ptr<SweepChunk> sweep_chunk)
  : WGSContext(lbs_solver, groupset, set_source_function, lhs_scope, rhs_scope, log_info),
    sweep_chunk_(std::move(sweep_chunk)),
    sweep_scheduler_(lbs_solver.SweepType() == "AAH" ? SchedulingAlgorithm::DEPTH_OF_GRAPH
                                                     : SchedulingAlgorithm::FIRST_IN_FIRST_OUT,
                     *groupset.angle_agg_,
                     *sweep_chunk_),
    lbs_ss_solver_(lbs_solver)
{
}

void
SweepWGSContext::PreSetupCallback()
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::PreSetupCallback");

  if (log_info_)
  {
    std::string method_name;
    switch (groupset_.iterative_method_)
    {
      case IterativeMethod::KRYLOV_RICHARDSON:
        method_name = "KRYLOV_RICHARDSON";
        break;
      case IterativeMethod::KRYLOV_GMRES:
        method_name = "KRYLOV_GMRES";
        break;
      case IterativeMethod::KRYLOV_BICGSTAB:
        method_name = "KRYLOV_BICGSTAB";
        break;
      default:
        method_name = "KRYLOV_GMRES";
    }
    log.Log() << "\n\n"
              << "********** Solving groupset " << groupset_.id_ << " with " << method_name
              << ".\n\n"
              << "Quadrature number of angles: " << groupset_.quadrature_->abscissae_.size() << "\n"
              << "Groups " << groupset_.groups_.front().id_ << " " << groupset_.groups_.back().id_
              << "\n\n";
  }
}

void
SweepWGSContext::SetPreconditioner(KSP& solver)
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::SetPreconditioner");

  auto& ksp = solver;

  PC pc;
  KSPGetPC(ksp, &pc);

  if (groupset_.apply_wgdsa_ or groupset_.apply_tgdsa_)
  {
    PCSetType(pc, PCSHELL);
    PCShellSetApply(pc, (PCShellPtr)WGDSA_TGDSA_PreConditionerMult);
    PCShellSetContext(pc, &(*this));
  }

  KSPSetPCSide(ksp, PC_LEFT);
  KSPSetUp(ksp);
}

std::pair<int64_t, int64_t>
SweepWGSContext::SystemSize()
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::SystemSize");

  const size_t local_node_count = lbs_solver_.LocalNodeCount();
  const size_t globl_node_count = lbs_solver_.GlobalNodeCount();
  const size_t num_moments = lbs_solver_.NumMoments();

  const size_t groupset_numgrps = groupset_.groups_.size();
  const auto num_delayed_psi_info = groupset_.angle_agg_->GetNumDelayedAngularDOFs();
  const size_t local_size =
    local_node_count * num_moments * groupset_numgrps + num_delayed_psi_info.first;
  const size_t globl_size =
    globl_node_count * num_moments * groupset_numgrps + num_delayed_psi_info.second;
  const size_t num_angles = groupset_.quadrature_->abscissae_.size();
  const size_t num_psi_global = globl_node_count * num_angles * groupset_.groups_.size();
  const size_t num_delayed_psi_globl = num_delayed_psi_info.second;

  if (log_info_)
  {
    log.Log() << "Total number of angular unknowns: " << num_psi_global << "\n"
              << "Number of lagged angular unknowns: " << num_delayed_psi_globl << "("
              << std::setprecision(2)
              << static_cast<double>(num_delayed_psi_globl) * 100 /
                   static_cast<double>(num_psi_global)
              << "%)";
  }

  return {static_cast<int64_t>(local_size), static_cast<int64_t>(globl_size)};
}

void
SweepWGSContext::ApplyInverseTransportOperator(SourceFlags scope)
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::ApplyInverseTransportOperator");

  ++counter_applications_of_inv_op_;
  const bool use_bndry_source_flag =
    (scope & APPLY_FIXED_SOURCES) and (not lbs_solver_.Options().use_src_moments);

  sweep_scheduler_.SetBoundarySourceActiveFlag(use_bndry_source_flag);

  if (scope & ZERO_INCOMING_DELAYED_PSI)
    sweep_scheduler_.ZeroIncomingDelayedPsi();

  // Sweep
  sweep_scheduler_.ZeroOutputFluxDataStructures();
  std::chrono::high_resolution_clock::time_point sweep_start =
    std::chrono::high_resolution_clock::now();
  sweep_scheduler_.Sweep();
  std::chrono::high_resolution_clock::time_point sweep_end =
    std::chrono::high_resolution_clock::now();
  double sweep_time =
    (std::chrono::duration_cast<std::chrono::nanoseconds>(sweep_end - sweep_start).count()) /
    1.0e+9;
  sweep_times_.push_back(sweep_time);
}

void
SweepWGSContext::PostSolveCallback()
{
  CALI_CXX_MARK_SCOPE("SweepWGSContext::PostSolveCallback");

  // Perform final sweep with converged phi and delayed psi dofs
  if (groupset_.iterative_method_ != IterativeMethod::KRYLOV_RICHARDSON)
  {
    lbs_ss_solver_.ZeroOutflowBalanceVars(groupset_);

    const auto scope = lhs_src_scope_ | rhs_src_scope_;

    set_source_function_(groupset_,
                         lbs_solver_.QMomentsLocal(),
                         lbs_solver_.PhiOldLocal(),
                         lbs_solver_.DensitiesLocal(),
                         scope);
    sweep_scheduler_.SetDestinationPhi(lbs_solver_.PhiNewLocal());

    ApplyInverseTransportOperator(scope);

    lbs_solver_.GSScopedCopyPrimarySTLvectors(
      groupset_, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);
  }

  double tot_sweep_time = 0.0;
  for (auto time : sweep_times_)
    tot_sweep_time += time;
  double num_sweeps = static_cast<double>(sweep_times_.size());
  double avg_sweep_time = tot_sweep_time / num_sweeps;
  size_t num_angles = groupset_.quadrature_->abscissae_.size();
  size_t num_unknowns = lbs_solver_.GlobalNodeCount() * num_angles * groupset_.groups_.size();

  log.Log() << "\n       Average sweep time (s):        "
            << tot_sweep_time / static_cast<double>(sweep_times_.size())
            << "\n       Sweep Time/Unknown (ns):       "
            << avg_sweep_time * 1.0e9 * opensn::mpi_comm.size() / static_cast<double>(num_unknowns)
            << "\n       Number of unknowns per sweep:  " << num_unknowns << "\n\n";
}

} // namespace lbs
} // namespace opensn
