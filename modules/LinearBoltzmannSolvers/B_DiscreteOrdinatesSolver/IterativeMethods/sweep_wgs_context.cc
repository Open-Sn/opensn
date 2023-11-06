#include "modules/LinearBoltzmannSolvers/B_DiscreteOrdinatesSolver/IterativeMethods/sweep_wgs_context.h"

#include <petscksp.h>

#include "modules/LinearBoltzmannSolvers/B_DiscreteOrdinatesSolver/lbs_discrete_ordinates_solver.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Preconditioning/lbs_shell_operations.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include <iomanip>

#define sc_double static_cast<double>
#define PCShellPtr PetscErrorCode (*)(PC, Vec, Vec)

namespace lbs
{

/**General print out of information.*/
template <>
void
SweepWGSContext<Mat, Vec, KSP>::PreSetupCallback()
{
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
    Chi::log.Log() << "\n\n"
                   << "********** Solving groupset " << groupset_.id_ << " with " << method_name
                   << ".\n\n"
                   << "Quadrature number of angles: " << groupset_.quadrature_->abscissae_.size()
                   << "\n"
                   << "Groups " << groupset_.groups_.front().id_ << " "
                   << groupset_.groups_.back().id_ << "\n\n";
  }
}

/**Sets the preconditioner application function.*/
template <>
void
SweepWGSContext<Mat, Vec, KSP>::SetPreconditioner(KSP& solver)
{
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

/**For sweeping we add lagged angular fluxes to the size of the vectors.*/
template <>
std::pair<int64_t, int64_t>
SweepWGSContext<Mat, Vec, KSP>::SystemSize()
{
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
    Chi::log.Log() << "Total number of angular unknowns: " << num_psi_global << "\n"
                   << "Number of lagged angular unknowns: " << num_delayed_psi_globl << "("
                   << std::setprecision(2)
                   << sc_double(num_delayed_psi_globl) * 100 / sc_double(num_psi_global) << "%)";
  }

  return {static_cast<int64_t>(local_size), static_cast<int64_t>(globl_size)};
}

/**With a right-hand side built. This routine applies the inverse
 * of the transport operator to this right-hand side.*/
template <>
void
SweepWGSContext<Mat, Vec, KSP>::ApplyInverseTransportOperator(int scope)
{
  ++counter_applications_of_inv_op_;
  const bool use_bndry_source_flag =
    (scope & APPLY_FIXED_SOURCES) and (not lbs_solver_.Options().use_src_moments);

  sweep_scheduler_.SetBoundarySourceActiveFlag(use_bndry_source_flag);

  if (scope & ZERO_INCOMING_DELAYED_PSI) sweep_scheduler_.ZeroIncomingDelayedPsi();

  // Sweep
  sweep_scheduler_.ZeroOutputFluxDataStructures();
  sweep_scheduler_.Sweep();
}

/**This method implements an additional sweep for two reasons:
 * The first is to compute balance parameters, and the second
 * is to allow for the calculation of proper angular fluxes. The
 * latter is needed because some krylov methods do not necessarily
 * provide the true angular flux at each iteration.*/
template <>
void
SweepWGSContext<Mat, Vec, KSP>::PostSolveCallback()
{
  //================================================== Perform final sweep
  //                                                   with converged phi and
  //                                                   delayed psi dofs
  if (groupset_.iterative_method_ != IterativeMethod::KRYLOV_RICHARDSON)
  {
    lbs_ss_solver_.ZeroOutflowBalanceVars(groupset_);

    const int scope = lhs_src_scope_ | rhs_src_scope_;

    set_source_function_(groupset_, lbs_solver_.QMomentsLocal(), lbs_solver_.PhiOldLocal(), scope);
    sweep_scheduler_.SetDestinationPhi(lbs_solver_.PhiNewLocal());

    ApplyInverseTransportOperator(scope);

    lbs_solver_.GSScopedCopyPrimarySTLvectors(
      groupset_, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);
  }

  //==================================================== Print solution info
  {
    double sweep_time = sweep_scheduler_.GetAverageSweepTime();
    double chunk_overhead_ratio = 1.0 - sweep_scheduler_.GetAngleSetTimings()[2];
    double source_time = Chi::log.ProcessEvent(lbs_solver_.GetSourceEventTag(),
                                               chi::ChiLog::EventOperation::AVERAGE_DURATION);
    size_t num_angles = groupset_.quadrature_->abscissae_.size();
    size_t num_unknowns = lbs_solver_.GlobalNodeCount() * num_angles * groupset_.groups_.size();

    if (log_info_)
    {
      Chi::log.Log() << "\n\n";
      Chi::log.Log() << "        Set Src Time/sweep (s):        " << source_time;
      Chi::log.Log() << "        Average sweep time (s):        " << sweep_time;
      Chi::log.Log() << "        Chunk-Overhead-Ratio  :        " << chunk_overhead_ratio;
      Chi::log.Log() << "        Sweep Time/Unknown (ns):       "
                     << sweep_time * 1.0e9 * Chi::mpi.process_count /
                          static_cast<double>(num_unknowns);
      Chi::log.Log() << "        Number of unknowns per sweep:  " << num_unknowns;
      Chi::log.Log() << "\n\n";

      std::string sweep_log_file_name = std::string("GS_") + std::to_string(groupset_.id_) +
                                        std::string("_SweepLog_") +
                                        std::to_string(Chi::mpi.location_id) + std::string(".log");
      groupset_.PrintSweepInfoFile(sweep_scheduler_.SweepEventTag(), sweep_log_file_name);
    }
  }
}

} // namespace lbs
