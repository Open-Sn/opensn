// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"
#include "framework/math/linear_solver/petsc_linear_solver.h"
#include <memory>
#include <vector>
#include <functional>

namespace opensn
{

/// Linear Solver specialization for Within GroupSet (WGS) solves.
class WGSLinearSolver : public PETScLinearSolver
{
public:
  /**
   * Constructor.
   * \param gs_context_ptr Context Pointer to abstract context.
   */
  explicit WGSLinearSolver(const std::shared_ptr<WGSContext>& gs_context_ptr);

  ~WGSLinearSolver() override;

protected:
  void PreSetupCallback() override;
  void SetConvergenceTest() override;
  void SetSystemSize() override;
  void SetSystem() override;
  void SetPreconditioner() override;
  void PostSetupCallback() override;
  void PreSolveCallback() override;
  void SetRHS() override;
  void SetInitialGuess() override;
  void PostSolveCallback() override;

  std::vector<double> saved_q_moments_local_;
};

} // namespace opensn
