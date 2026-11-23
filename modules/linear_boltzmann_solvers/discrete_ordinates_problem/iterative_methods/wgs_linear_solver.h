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

  /**
   * Builds the system if needed and solves it. Rebuilds the system first when x_ and b_ view flux
   * storage that has been reallocated since they were built.
   */
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

  /// Work vector for the preconditioned right-hand side; a view of the context's source scratch.
  Vec rhs_preconditioned_work_;
  /**
   * True when x_ and b_ view the groupset's owned blocks of phi_old and phi_new, extended into the
   * vectors' tails by the local delayed angular unknowns. That is possible when the groupset has no
   * delayed angular unknowns on any rank, or when its block is the last one. Otherwise x_ and b_
   * own their storage. The views are valid for the generations in phi_old_generation_ and
   * phi_new_generation_.
   */
  bool aliases_flux_ = false;
  std::uint64_t phi_old_generation_ = 0;
  std::uint64_t phi_new_generation_ = 0;
};

} // namespace opensn
