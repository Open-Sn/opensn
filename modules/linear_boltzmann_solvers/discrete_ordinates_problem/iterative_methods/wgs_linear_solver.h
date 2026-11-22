// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "framework/math/linear_solver/petsc_linear_system_solver.h"
#include <cstdint>
#include <memory>

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

  void Solve() override;

  /// Destroys the KSP, the operator, the vectors, and the source scratch of the context.
  void ReleaseWorkspace() override;

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

  Vec rhs_preconditioned_work_;
  /**
   * True when x_ is a view of the groupset's owned phi_old block, which is possible when the
   * groupset has no delayed angular unknowns. The view is valid for the generation of phi_old in
   * x_alias_generation_.
   */
  bool x_aliases_phi_old_ = false;
  std::uint64_t x_alias_generation_ = 0;
};

} // namespace opensn
