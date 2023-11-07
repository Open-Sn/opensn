#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Preconditioning/lbs_shell_operations.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/diffusion_mip_solver.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_context.h"

/**Applies TGDSA to the given input vector.*/
int
lbs::MIP_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output)
{
  void* context;
  PCShellGetContext(pc, &context);

  auto gs_context_ptr = (lbs::WGSContext<Mat, Vec, KSP>*)(context);

  // Shorten some names
  lbs::LBSSolver& solver = gs_context_ptr->lbs_solver_;
  LBSGroupset& groupset = gs_context_ptr->groupset_;

  // Copy PETSc vector to STL
  auto& phi_delta = gs_context_ptr->lbs_solver_.PhiNewLocal();
  solver.SetPrimarySTLvectorFromGSPETScVec(groupset, phi_input, PhiSTLOption::PHI_NEW);

  // Apply TGDSA
  if (groupset.apply_tgdsa_)
  {
    std::vector<double> delta_phi_local;
    solver.AssembleTGDSADeltaPhiVector(groupset, phi_delta, delta_phi_local);
    groupset.tgdsa_solver_->Assemble_b(delta_phi_local);
    groupset.tgdsa_solver_->Solve(delta_phi_local);
    solver.DisAssembleTGDSADeltaPhiVector(groupset, delta_phi_local, phi_delta);
  }

  // Copy STL vector to PETSc Vec
  solver.SetGSPETScVecFromPrimarySTLvector(groupset, pc_output, PhiSTLOption::PHI_NEW);

  return 0;
}
