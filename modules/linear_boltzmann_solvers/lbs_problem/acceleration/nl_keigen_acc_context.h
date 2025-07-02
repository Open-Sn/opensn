// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/snes_k_residual_func_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/math/nonlinear_solver/nonlinear_solver_context.h"
#include <petscsnes.h>

namespace opensn
{

struct NLKEigenDiffContext : public NonLinearSolverContext
{
  DiffusionMIPSolver& diff_solver;
  LBSProblem& lbs_problem;
  int verbosity_level;
  KResidualFunctionContext kresid_func_context;

  size_t diff_num_local_dofs;

  std::vector<double> phi_l;
  std::vector<double> phi_lph_i;
  std::vector<double> phi_lph_ip1;
  std::vector<double> Sf;
  double k_l = 1.0;

  explicit NLKEigenDiffContext(DiffusionMIPSolver& diff_solver,
                               LBSProblem& lbs_problem,
                               int verbosity_level)
    : diff_solver(diff_solver),
      lbs_problem(lbs_problem),
      verbosity_level(verbosity_level),
      kresid_func_context({diff_solver.GetName(), 1.0}),
      diff_num_local_dofs(diff_solver.GetNumPhiIterativeUnknowns().first)
  {
    phi_l.assign(diff_num_local_dofs, 0.0);
    phi_lph_i.assign(diff_num_local_dofs, 0.0);
    phi_lph_ip1.assign(diff_num_local_dofs, 0.0);
    Sf.assign(diff_num_local_dofs, 0.0);
  }

  std::vector<double> PhiVecToSTLVec(Vec phi) const
  {
    std::vector<double> output(diff_num_local_dofs, 0.0);

    const double* phi_raw;
    VecGetArrayRead(phi, &phi_raw);
    for (size_t i = 0; i < diff_num_local_dofs; ++i)
      output[i] = phi_raw[i];
    VecRestoreArrayRead(phi, &phi_raw);

    return output;
  }

  void STLVecToPhiVec(const std::vector<double>& input, Vec phi) const
  {
    double* phi_raw;
    VecGetArray(phi, &phi_raw);
    for (size_t i = 0; i < diff_num_local_dofs; ++i)
      phi_raw[i] = input[i];
    VecRestoreArray(phi, &phi_raw);
  }

  ~NLKEigenDiffContext() override = default;
};

} // namespace opensn
