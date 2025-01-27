// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/diffusion_dfem_solver/iterative_methods/mip_wgs_context.h"
#include "modules/linear_boltzmann_solvers/diffusion_dfem_solver/lbs_mip_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/preconditioning/lbs_shell_operations.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion_mip_solver.h"
#include "framework/runtime.h"
#include <petscksp.h>
#include <iomanip>

namespace opensn
{

using PCShellPtr = PetscErrorCode (*)(PC, Vec, Vec);

MIPWGSContext::MIPWGSContext(DiffusionDFEMSolver& solver,
                             LBSGroupset& groupset,
                             const SetSourceFunction& set_source_function,
                             SourceFlags lhs_scope,
                             SourceFlags rhs_scope,
                             bool log_info)
  : WGSContext(solver, groupset, set_source_function, lhs_scope, rhs_scope, log_info)
{
}

void
MIPWGSContext::PreSetupCallback()
{
  if (log_info)
    log.Log() << "\n\n"
              << "********** Solving groupset " << groupset.id << " with "
              << LinearSolver::IterativeMethodName(groupset.iterative_method) << ".\n\n";
}

void
MIPWGSContext::SetPreconditioner(KSP& solver)
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
MIPWGSContext::GetSystemSize()
{
  const size_t local_node_count = lbs_solver.GetLocalNodeCount();
  const size_t globl_node_count = lbs_solver.GetGlobalNodeCount();

  const size_t groupset_numgrps = groupset.groups.size();
  const size_t local_size = local_node_count * groupset_numgrps;
  const size_t globl_size = globl_node_count * groupset_numgrps;

  return {static_cast<int64_t>(local_size), static_cast<int64_t>(globl_size)};
}

void
MIPWGSContext::ApplyInverseTransportOperator(SourceFlags scope)
{
  ++counter_applications_of_inv_op;
  auto& mip_solver = *dynamic_cast<DiffusionDFEMSolver&>(lbs_solver).gs_mip_solvers[groupset.id];

  lbs_solver.GetPhiNewLocal() = lbs_solver.GetQMomentsLocal();

  Vec work_vector;
  VecDuplicate(mip_solver.GetRHS(), &work_vector);

  LBSVecOps::SetGSPETScVecFromPrimarySTLvector(
    lbs_solver, groupset, work_vector, PhiSTLOption::PHI_NEW);

  mip_solver.Assemble_b(work_vector);
  mip_solver.Solve(work_vector);

  LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(
    lbs_solver, groupset, work_vector, PhiSTLOption::PHI_NEW);

  VecDestroy(&work_vector);
}

void
MIPWGSContext::PostSolveCallback()
{
  LBSVecOps::GSScopedCopyPrimarySTLvectors(
    lbs_solver, groupset, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);
}

} // namespace opensn
