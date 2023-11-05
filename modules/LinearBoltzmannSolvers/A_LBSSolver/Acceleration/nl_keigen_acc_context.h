#pragma once

#include "opensn/framework/math/NonLinearSolver/NonLinearSolverContext.h"

#include "opensn/modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/diffusion_mip_solver.h"
#include "opensn/modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/snes_k_residual_func_context.h"
#include "opensn/modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include <petscsnes.h>

namespace lbs::acceleration
{

struct NLKEigenDiffContext : public chi_math::NonLinearSolverContext<Vec, SNES>
{
  DiffusionMIPSolver& diff_solver_;
  LBSSolver& lbs_solver_;
  int verbosity_level_;
  KResidualFunctionContext kresid_func_context_;

  size_t diff_num_local_dofs_;

  VecDbl phi_l_;
  VecDbl phi_lph_i_;
  VecDbl phi_lph_ip1_;
  VecDbl Sf_;
  double k_l = 1.0;

  explicit NLKEigenDiffContext(DiffusionMIPSolver& diff_solver,
                               LBSSolver& lbs_solver,
                               int verbosity_level)
    : diff_solver_(diff_solver),
      lbs_solver_(lbs_solver),
      verbosity_level_(verbosity_level),
      kresid_func_context_({diff_solver.TextName(), 1.0}),
      diff_num_local_dofs_(diff_solver_.GetNumPhiIterativeUnknowns().first)
  {
  }

  VecDbl PhiVecToSTLVec(Vec phi) const
  {
    VecDbl output(diff_num_local_dofs_, 0.0);

    const double* phi_raw;
    VecGetArrayRead(phi, &phi_raw);
    for (size_t i = 0; i < diff_num_local_dofs_; ++i)
      output[i] = phi_raw[i];
    VecRestoreArrayRead(phi, &phi_raw);

    return output;
  }

  void STLVecToPhiVec(const VecDbl& input, Vec phi) const
  {
    double* phi_raw;
    VecGetArray(phi, &phi_raw);
    for (size_t i = 0; i < diff_num_local_dofs_; ++i)
      phi_raw[i] = input[i];
    VecRestoreArray(phi, &phi_raw);
  }

  virtual ~NLKEigenDiffContext() override = default;
};

} // namespace lbs::acceleration
