// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_system_solver.h"
#include "framework/math/linear_solver/linear_solver_context.h"
#include <petscsystypes.h>
#include <string>
#include <utility>
#include <memory>
#include <petscksp.h>

namespace opensn
{

class PETScLinearSolver : public LinearSystemSolver
{
public:
  struct ToleranceOptions
  {
    double residual_relative = 1.0e-50;
    double residual_absolute = 1.0e-6;
    double residual_divergence = 1.0e6;
    PetscInt maximum_iterations = 100;
    PetscInt gmres_restart_interval = 100;
    double gmres_breakdown_tolerance = 1.0e6;
  } tolerance_options;

  PETScLinearSolver(IterativeMethod method, std::shared_ptr<LinearSystemContext> context_ptr);

  ~PETScLinearSolver() override;

  ToleranceOptions& GetToleranceOptions() { return tolerance_options; }

  void ApplyToleranceOptions();

  /// Sets a flag to suppress the KSPSolve() method from being called.
  void SetKSPSolveSuppressionFlag(bool flag) { suppress_kspsolve_ = flag; }

  bool GetKSPSolveSuppressionFlag() const { return suppress_kspsolve_; }

  /// Set up the linaer solver
  void Setup() override;

  /// Solve the system
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
  virtual void PostSetupCallback();
  virtual void PreSolveCallback();
  virtual void SetInitialGuess() = 0;
  virtual void SetRHS() = 0;
  virtual void PostSolveCallback();

  Mat A_;
  Vec b_;
  Vec x_;
  KSP ksp_;
  int64_t num_local_dofs_;
  int64_t num_global_dofs_;

private:
  bool system_set_;
  bool suppress_kspsolve_;
  std::string PETScIterativeMethodName();

protected:
  static PetscErrorCode LinearSolverMatrixAction(Mat matrix, Vec vector, Vec action);
};

} // namespace opensn
