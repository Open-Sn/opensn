// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/preconditioning/lbs_shell_operations.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/tgdsa.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "caliper/cali.h"
#include <cstdint>
#include <span>
#include <vector>

namespace opensn
{

PetscErrorCode
MIP_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output)
{
  void* context = nullptr;
  PetscErrorCode ierr = PCShellGetContext(pc, static_cast<void*>(&context));
  if (ierr != PETSC_SUCCESS)
    return ierr;

  auto* gs_context_ptr = static_cast<WGSContext*>(context);

  // Shorten some names
  DiscreteOrdinatesProblem& solver = gs_context_ptr->do_problem;
  LBSGroupset& groupset = gs_context_ptr->groupset;

  // TGDSA corrects the flux part of the output in place.
  const auto block_size = solver.GetPhiNewLocal().GetBlockSize(solver.GetFluxBlock(groupset));
  const auto last_flux_index = static_cast<int64_t>(block_size) - 1;
  if (phi_input != pc_output)
    OpenSnPETScCall(VecCopy(phi_input, pc_output));
  PetscScalar* output = nullptr;
  OpenSnPETScCall(VecGetArray(pc_output, &output));
  const std::span<double> phi_delta(output, block_size);
  if (groupset.angle_agg)
  {
    auto index = last_flux_index;
    groupset.angle_agg->SetNewDelayedAngularDOFsFromArray(index, output);
  }

  // Apply TGDSA
  if (groupset.apply_tgdsa)
  {
    CALI_CXX_MARK_SCOPE("Acceleration/TGDSA");

    // DSA dsa(solver);
    std::vector<double> delta_phi_local;
    TGDSA::AssembleDeltaPhiVector(solver, groupset, phi_delta, delta_phi_local);
    groupset.tgdsa_solver->Assemble_b(delta_phi_local);
    groupset.tgdsa_solver->Solve(delta_phi_local);
    TGDSA::DisassembleDeltaPhiVector(solver, groupset, delta_phi_local, phi_delta);
  }

  if (groupset.angle_agg)
  {
    auto index = last_flux_index;
    groupset.angle_agg->AppendNewDelayedAngularDOFsToArray(index, output);
  }
  OpenSnPETScCall(VecRestoreArray(pc_output, &output));

  return PETSC_SUCCESS;
}

} // namespace opensn
