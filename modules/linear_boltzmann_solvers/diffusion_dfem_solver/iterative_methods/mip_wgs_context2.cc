// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/diffusion_dfem_solver/iterative_methods/mip_wgs_context2.h"
#include "modules/linear_boltzmann_solvers/diffusion_dfem_solver/lbs_mip_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/preconditioning/lbs_shell_operations.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion_mip_solver.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <petscksp.h>
#include <iomanip>

namespace opensn
{

using PCShellPtr = PetscErrorCode (*)(PC, Vec, Vec);

MIPWGSContext2::MIPWGSContext2(DiffusionDFEMSolver& lbs_mip_ss_solver,
                               LBSGroupset& groupset,
                               const SetSourceFunction& set_source_function,
                               SourceFlags lhs_scope,
                               SourceFlags rhs_scope,
                               bool log_info)
  : WGSContext(lbs_mip_ss_solver, groupset, set_source_function, lhs_scope, rhs_scope, log_info),
    lbs_mip_ss_solver(lbs_mip_ss_solver)
{
}

void
MIPWGSContext2::PreSetupCallback()
{
  if (log_info)
  {
    std::string method_name;
    switch (groupset.iterative_method)
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
              << "********** Solving groupset " << groupset.id << " with " << method_name
              << ".\n\n";
  }
}

void
MIPWGSContext2::SetPreconditioner(KSP& solver)
{
  auto& ksp = solver;

  PC pc;
  KSPGetPC(ksp, &pc);

  if (groupset.apply_tgdsa)
  {
    PCSetType(pc, PCSHELL);
    PCShellSetApply(pc, (PCShellPtr)MIP_TGDSA_PreConditionerMult);
    PCShellSetContext(pc, &(*this));
  }

  KSPSetPCSide(ksp, PC_LEFT);
  KSPSetUp(ksp);
}

std::pair<int64_t, int64_t>
MIPWGSContext2::SystemSize()
{
  const size_t local_node_count = lbs_solver.LocalNodeCount();
  const size_t globl_node_count = lbs_solver.GlobalNodeCount();

  const size_t groupset_numgrps = groupset.groups.size();
  const size_t local_size = local_node_count * groupset_numgrps;
  const size_t globl_size = globl_node_count * groupset_numgrps;

  return {static_cast<int64_t>(local_size), static_cast<int64_t>(globl_size)};
}

void
MIPWGSContext2::ApplyInverseTransportOperator(SourceFlags scope)
{
  ++counter_applications_of_inv_op;
  auto& mip_solver = *lbs_mip_ss_solver.gs_mip_solvers[groupset.id];

  lbs_solver.PhiNewLocal() = lbs_solver.QMomentsLocal();

  Vec work_vector;
  VecDuplicate(mip_solver.RHS(), &work_vector);

  lbs_solver.SetGSPETScVecFromPrimarySTLvector(groupset, work_vector, PhiSTLOption::PHI_NEW);

  mip_solver.Assemble_b(work_vector);
  mip_solver.Solve(work_vector);

  lbs_solver.SetPrimarySTLvectorFromGSPETScVec(groupset, work_vector, PhiSTLOption::PHI_NEW);

  VecDestroy(&work_vector);
}

void
MIPWGSContext2::PostSolveCallback()
{
  lbs_solver.GSScopedCopyPrimarySTLvectors(groupset, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);
}

} // namespace opensn
