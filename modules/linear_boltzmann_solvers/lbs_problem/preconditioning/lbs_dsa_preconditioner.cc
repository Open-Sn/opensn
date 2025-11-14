// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/preconditioning/lbs_shell_operations.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/wgdsa.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/tgdsa.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/diffusion/diffusion_mip_solver.h"

namespace opensn
{

PetscErrorCode
WGDSA_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output)
{
  void* context = nullptr;
  PCShellGetContext(pc, static_cast<void*>(&context));

  auto* gs_context_ptr = (WGSContext*)(context);

  // Shorten some names
  DiscreteOrdinatesProblem& do_problem = gs_context_ptr->do_problem;
  LBSGroupset& groupset = gs_context_ptr->groupset;

  // Copy PETSc vector to STL
  auto& phi_new_local = gs_context_ptr->do_problem.GetPhiNewLocal();
  LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(
    do_problem, groupset, phi_input, PhiSTLOption::PHI_NEW);

  // Apply WGDSA
  if (groupset.apply_wgdsa)
  {
    std::vector<double> delta_phi_local;
    WGDSA::AssembleDeltaPhiVector(do_problem, groupset, phi_new_local, delta_phi_local);

    groupset.wgdsa_solver->Assemble_b(delta_phi_local);
    groupset.wgdsa_solver->Solve(delta_phi_local);

    WGDSA::DisassembleDeltaPhiVector(do_problem, groupset, delta_phi_local, phi_new_local);
  }
  // Apply TGDSA
  if (groupset.apply_tgdsa)
  {
    std::vector<double> delta_phi_local;
    TGDSA::AssembleDeltaPhiVector(do_problem, groupset, phi_new_local, delta_phi_local);

    groupset.tgdsa_solver->Assemble_b(delta_phi_local);
    groupset.tgdsa_solver->Solve(delta_phi_local);

    TGDSA::DisassembleDeltaPhiVector(do_problem, groupset, delta_phi_local, phi_new_local);
  }

  // Copy STL vector to PETSc Vec
  LBSVecOps::SetGSPETScVecFromPrimarySTLvector(
    do_problem, groupset, pc_output, PhiSTLOption::PHI_NEW);

  return 0;
}

PetscErrorCode
WGDSA_TGDSA_PreConditionerMult2(WGSContext& gs_context_ptr, Vec phi_input, Vec pc_output)
{
  // Shorten some names
  DiscreteOrdinatesProblem& do_problem = gs_context_ptr.do_problem;
  LBSGroupset& groupset = gs_context_ptr.groupset;

  // Copy PETSc vector to STL
  auto& phi_new_local = gs_context_ptr.do_problem.GetPhiNewLocal();
  LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(
    do_problem, groupset, phi_input, PhiSTLOption::PHI_NEW);

  // Apply WGDSA
  if (groupset.apply_wgdsa)
  {
    std::vector<double> delta_phi_local;
    WGDSA::AssembleDeltaPhiVector(do_problem, groupset, phi_new_local, delta_phi_local);

    groupset.wgdsa_solver->Assemble_b(delta_phi_local);
    groupset.wgdsa_solver->Solve(delta_phi_local);

    WGDSA::DisassembleDeltaPhiVector(do_problem, groupset, delta_phi_local, phi_new_local);
  }
  // Apply TGDSA
  if (groupset.apply_tgdsa)
  {
    std::vector<double> delta_phi_local;
    TGDSA::AssembleDeltaPhiVector(do_problem, groupset, phi_new_local, delta_phi_local);

    groupset.tgdsa_solver->Assemble_b(delta_phi_local);
    groupset.tgdsa_solver->Solve(delta_phi_local);

    TGDSA::DisassembleDeltaPhiVector(do_problem, groupset, delta_phi_local, phi_new_local);
  }

  // Copy STL vector to PETSc Vec
  LBSVecOps::SetGSPETScVecFromPrimarySTLvector(
    do_problem, groupset, pc_output, PhiSTLOption::PHI_NEW);

  return 0;
}

} // namespace opensn
