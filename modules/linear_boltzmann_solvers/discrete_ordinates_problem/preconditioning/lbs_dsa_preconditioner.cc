// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/preconditioning/lbs_shell_operations.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/wgdsa.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/tgdsa.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

namespace opensn
{

namespace
{

/**
 * Applies WGDSA/TGDSA to the groupset vector `phi_input` (flux moments followed by the lagged
 * angular unknowns) and writes the result to `pc_output`, which may be `phi_input`. The flux part
 * is accelerated in place in `pc_output`. The scalar-flux correction is also applied, as an
 * isotropic angular correction, to the lagged boundary angular fluxes.
 */
void
ApplyDSAPreconditioner(WGSContext& gs_context, Vec phi_input, Vec pc_output)
{
  DiscreteOrdinatesProblem& do_problem = gs_context.do_problem;
  LBSGroupset& groupset = gs_context.groupset;

  const auto block_size =
    do_problem.GetPhiNewLocal().GetBlockSize(do_problem.GetFluxBlock(groupset));
  // Index of the last flux entry; the delayed angular unknowns follow it.
  const auto last_flux_index = static_cast<int64_t>(block_size) - 1;
  const bool in_place = phi_input == pc_output;

  PetscScalar* output = nullptr;
  const PetscScalar* input = nullptr;
  OpenSnPETScCall(VecGetArray(pc_output, &output));
  if (in_place)
    input = output;
  else
    OpenSnPETScCall(VecGetArrayRead(phi_input, &input));

  const std::span<double> phi(output, block_size);
  if (not in_place)
    std::copy_n(input, block_size, phi.begin());
  // The lagged angular fluxes of the input, which the WGDSA source and the correction use.
  if (groupset.angle_agg)
  {
    auto index = last_flux_index;
    groupset.angle_agg->SetNewDelayedAngularDOFsFromArray(index, input);
  }

  const bool correct_delayed_psi =
    groupset.angle_agg and groupset.angle_agg->GetNumDelayedAngularDOFs().second > 0;
  // Change of the flux by the acceleration. In place, the input is overwritten, so keep it.
  std::vector<double> delta_phi;
  if (correct_delayed_psi and in_place)
    delta_phi.assign(phi.begin(), phi.end());

  // Apply WGDSA
  if (groupset.apply_wgdsa)
  {
    CALI_CXX_MARK_SCOPE("Acceleration/WGDSA");

    WGDSA::Apply(do_problem, groupset, phi, false, phi);
  }
  // Apply TGDSA
  if (groupset.apply_tgdsa)
  {
    CALI_CXX_MARK_SCOPE("Acceleration/TGDSA");

    std::vector<double> delta_phi_local;
    TGDSA::AssembleDeltaPhiVector(do_problem, groupset, phi, delta_phi_local);

    groupset.tgdsa_solver->Assemble_b(delta_phi_local);
    groupset.tgdsa_solver->Solve(delta_phi_local);

    TGDSA::DisassembleDeltaPhiVector(do_problem, groupset, delta_phi_local, phi);
  }

  if (correct_delayed_psi)
  {
    if (in_place)
      for (size_t i = 0; i < block_size; ++i)
        delta_phi[i] = phi[i] - delta_phi[i];
    else
    {
      delta_phi.resize(block_size);
      for (size_t i = 0; i < block_size; ++i)
        delta_phi[i] = phi[i] - input[i];
    }
    ApplyDSACorrectionToDelayedAngularFlux(do_problem, groupset, delta_phi);
  }

  if (groupset.angle_agg)
  {
    auto index = last_flux_index;
    groupset.angle_agg->AppendNewDelayedAngularDOFsToArray(index, output);
  }

  if (not in_place)
    OpenSnPETScCall(VecRestoreArrayRead(phi_input, &input));
  OpenSnPETScCall(VecRestoreArray(pc_output, &output));
}

} // namespace

PetscErrorCode
WGDSA_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output)
{
  void* context = nullptr;
  PetscErrorCode ierr = PCShellGetContext(pc, static_cast<void*>(&context));
  if (ierr != PETSC_SUCCESS)
    return ierr;

  ApplyDSAPreconditioner(*static_cast<WGSContext*>(context), phi_input, pc_output);
  return PETSC_SUCCESS;
}

PetscErrorCode
WGDSA_TGDSA_PreConditionerMult2(WGSContext& gs_context_ptr, Vec phi_input, Vec pc_output)
{
  ApplyDSAPreconditioner(gs_context_ptr, phi_input, pc_output);
  return PETSC_SUCCESS;
}

} // namespace opensn
