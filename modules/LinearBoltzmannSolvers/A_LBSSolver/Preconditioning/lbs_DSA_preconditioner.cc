#include "lbs_shell_operations.h"

#include "A_LBSSolver/lbs_solver.h"
#include "A_LBSSolver/Acceleration/diffusion_mip.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_context.h"

//###################################################################
/**Applies WGDSA or TGDSA to the given input vector.*/
int
lbs::WGDSA_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output)
{
  void* context;
  PCShellGetContext(pc, &context);

  auto gs_context_ptr = (lbs::WGSContext<Mat, Vec, KSP>*)(context);

  // Shorten some names
  lbs::LBSSolver& lbs_solver = gs_context_ptr->lbs_solver_;
  LBSGroupset& groupset = gs_context_ptr->groupset_;

  //============================================= Copy PETSc vector to STL
  auto& phi_new_local = gs_context_ptr->lbs_solver_.PhiNewLocal();
  lbs_solver.SetPrimarySTLvectorFromGSPETScVec(groupset, phi_input, PhiSTLOption::PHI_NEW);

  //============================================= Apply WGDSA
  if (groupset.apply_wgdsa_)
  {
    std::vector<double> delta_phi_local;
    lbs_solver.AssembleWGDSADeltaPhiVector(groupset,
                                           phi_new_local,    // From
                                           delta_phi_local); // To

    groupset.wgdsa_solver_->Assemble_b(delta_phi_local);
    groupset.wgdsa_solver_->Solve(delta_phi_local);

    lbs_solver.DisAssembleWGDSADeltaPhiVector(groupset,
                                              delta_phi_local, // From
                                              phi_new_local);  // To
  }
  //============================================= Apply TGDSA
  if (groupset.apply_tgdsa_)
  {
    std::vector<double> delta_phi_local;
    lbs_solver.AssembleTGDSADeltaPhiVector(groupset,
                                           phi_new_local,    // From
                                           delta_phi_local); // To

    groupset.tgdsa_solver_->Assemble_b(delta_phi_local);
    groupset.tgdsa_solver_->Solve(delta_phi_local);

    lbs_solver.DisAssembleTGDSADeltaPhiVector(groupset,
                                              delta_phi_local, // From
                                              phi_new_local);  // To
  }

  //============================================= Copy STL vector to PETSc Vec
  lbs_solver.SetGSPETScVecFromPrimarySTLvector(groupset, pc_output, PhiSTLOption::PHI_NEW);

  return 0;
}

//###################################################################
/**Applies WGDSA or TGDSA to the given input vector.*/
int
lbs::WGDSA_TGDSA_PreConditionerMult2(lbs::WGSContext<Mat, Vec, KSP>& gs_context_ptr,
                                     Vec phi_input,
                                     Vec pc_output)
{
  // Shorten some names
  lbs::LBSSolver& lbs_solver = gs_context_ptr.lbs_solver_;
  LBSGroupset& groupset = gs_context_ptr.groupset_;

  //============================================= Copy PETSc vector to STL
  auto& phi_new_local = gs_context_ptr.lbs_solver_.PhiNewLocal();
  lbs_solver.SetPrimarySTLvectorFromGSPETScVec(groupset, phi_input, PhiSTLOption::PHI_NEW);

  //============================================= Apply WGDSA
  if (groupset.apply_wgdsa_)
  {
    std::vector<double> delta_phi_local;
    lbs_solver.AssembleWGDSADeltaPhiVector(groupset,
                                           phi_new_local,    // From
                                           delta_phi_local); // To

    groupset.wgdsa_solver_->Assemble_b(delta_phi_local);
    groupset.wgdsa_solver_->Solve(delta_phi_local);

    lbs_solver.DisAssembleWGDSADeltaPhiVector(groupset,
                                              delta_phi_local, // From
                                              phi_new_local);  // To
  }
  //============================================= Apply TGDSA
  if (groupset.apply_tgdsa_)
  {
    std::vector<double> delta_phi_local;
    lbs_solver.AssembleTGDSADeltaPhiVector(groupset,
                                           phi_new_local,    // From
                                           delta_phi_local); // To

    groupset.tgdsa_solver_->Assemble_b(delta_phi_local);
    groupset.tgdsa_solver_->Solve(delta_phi_local);

    lbs_solver.DisAssembleTGDSADeltaPhiVector(groupset,
                                              delta_phi_local, // From
                                              phi_new_local);  // To
  }

  //============================================= Copy STL vector to PETSc Vec
  lbs_solver.SetGSPETScVecFromPrimarySTLvector(groupset, pc_output, PhiSTLOption::PHI_NEW);

  return 0;
}