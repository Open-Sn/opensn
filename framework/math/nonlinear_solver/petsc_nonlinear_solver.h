// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/nonlinear_solver/nonlinear_solver.h"
#include "framework/math/nonlinear_solver/nonlinear_solver_context.h"
#include "framework/math/nonlinear_solver/petsc_nonlinear_solver_options.h"
#include <petscsnes.h>
#include <string>
#include <memory>
#include <utility>

namespace opensn
{

/// Implementation of a PETSc non-linear solver.
class PETScNonLinearSolver : public NonLinearSolver
{
public:
  explicit PETScNonLinearSolver(
    std::shared_ptr<NonLinearSolverContext> context_ptr,
    const InputParameters& params = PETScNonLinearSolverOptions::GetInputParameters());

  ~PETScNonLinearSolver() override;

  PETScNonLinearSolverOptions& GetToleranceOptions() { return options_; }

  void ApplyToleranceOptions();

  bool IsConverged() const { return converged_; }

  std::string GetConvergedReasonString() const;

  void Setup() override;

  void Solve() override;

protected:
  bool IsSystemSet() const { return system_set_; }
  virtual void PreSetupCallback();
  virtual void SetOptions();
  virtual void SetSolverContext();
  virtual void SetConvergenceTest();
  virtual void SetMonitor();
  virtual void SetPreconditioner();
  virtual void SetSystemSize() = 0;
  virtual void SetSystem() = 0;
  virtual void SetFunction() = 0;
  virtual void SetJacobian() = 0;
  virtual void PostSetupCallback();
  virtual void PreSolveCallback();
  virtual void SetInitialGuess() = 0;
  virtual void PostSolveCallback();

  Mat J_;
  Mat P_;
  Vec r_;
  Vec x_;
  SNES nl_solver_;

  int64_t num_local_dofs_;
  int64_t num_global_dofs_;

  PETScNonLinearSolverOptions options_;

private:
  bool system_set_;
  bool converged_;
  std::string converged_reason_string_;
};

} // namespace opensn
