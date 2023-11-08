#include "modules/LinearBoltzmannSolvers/B_DiffusionDFEMSolver/IterativeMethods/mip_wgs_context2.h"

#include <petscksp.h>

#include "modules/LinearBoltzmannSolvers/B_DiffusionDFEMSolver/lbsMIP_solver.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Preconditioning/lbs_shell_operations.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/diffusion_mip_solver.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include <iomanip>

#define sc_double static_cast<double>
#define PCShellPtr PetscErrorCode (*)(PC, Vec, Vec)

namespace lbs
{

template <>
void
MIPWGSContext2<Mat, Vec, KSP>::PreSetupCallback()
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
                   << ".\n\n";
  }
}

template <>
void
MIPWGSContext2<Mat, Vec, KSP>::SetPreconditioner(KSP& solver)
{
  auto& ksp = solver;

  PC pc;
  KSPGetPC(ksp, &pc);

  if (groupset_.apply_tgdsa_)
  {
    PCSetType(pc, PCSHELL);
    PCShellSetApply(pc, (PCShellPtr)MIP_TGDSA_PreConditionerMult);
    PCShellSetContext(pc, &(*this));
  }

  KSPSetPCSide(ksp, PC_LEFT);
  KSPSetUp(ksp);
}

template <>
std::pair<int64_t, int64_t>
MIPWGSContext2<Mat, Vec, KSP>::SystemSize()
{
  const size_t local_node_count = lbs_solver_.LocalNodeCount();
  const size_t globl_node_count = lbs_solver_.GlobalNodeCount();

  const size_t groupset_numgrps = groupset_.groups_.size();
  const size_t local_size = local_node_count * groupset_numgrps;
  const size_t globl_size = globl_node_count * groupset_numgrps;

  return {static_cast<int64_t>(local_size), static_cast<int64_t>(globl_size)};
}

template <>
void
MIPWGSContext2<Mat, Vec, KSP>::ApplyInverseTransportOperator(int scope)
{
  ++counter_applications_of_inv_op_;
  auto& mip_solver = *lbs_mip_ss_solver_.gs_mip_solvers_[groupset_.id_];

  lbs_solver_.PhiNewLocal() = lbs_solver_.QMomentsLocal();

  Vec work_vector;
  VecDuplicate(mip_solver.RHS(), &work_vector);

  lbs_solver_.SetGSPETScVecFromPrimarySTLvector(groupset_, work_vector, PhiSTLOption::PHI_NEW);

  mip_solver.Assemble_b(work_vector);
  mip_solver.Solve(work_vector);

  lbs_solver_.SetPrimarySTLvectorFromGSPETScVec(groupset_, work_vector, PhiSTLOption::PHI_NEW);

  VecDestroy(&work_vector);
}

template <>
void
MIPWGSContext2<Mat, Vec, KSP>::PostSolveCallback()
{
  lbs_solver_.GSScopedCopyPrimarySTLvectors(
    groupset_, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);
}

} // namespace lbs
