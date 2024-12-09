// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver.h"
#include "framework/math/linear_solver/linear_solver_context.h"
#include <string>
#include <utility>
#include <memory>
#include <petscksp.h>

namespace opensn
{

class PETScLinearSolver : public LinearSolver
{
public:
  struct ToleranceOptions
  {
    double residual_relative = 1.0e-50;
    double residual_absolute = 1.0e-6;
    double residual_divergence = 1.0e6;
    int maximum_iterations = 100;
    int gmres_restart_interval = 100;
    double gmres_breakdown_tolerance = 1.0e6;
  } tolerance_options;

  PETScLinearSolver(IterativeMethod iterative_method,
                    std::shared_ptr<LinearSolverContext> context_ptr);

  virtual ~PETScLinearSolver();

  ToleranceOptions& GetToleranceOptions() { return tolerance_options; }

  void ApplyToleranceOptions();

  /// Sets a flag to suppress the KSPSolve() method from being called.
  void SetKSPSolveSuppressionFlag(bool flag) { suppress_kspsolve_ = flag; }

  bool GetKSPSolveSuppressionFlag() const { return suppress_kspsolve_; }

  /// Set up the linaer solver
  virtual void Setup();

  /// Solve the system
  virtual void Solve();

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

protected:
  static int LinearSolverMatrixAction(Mat matrix, Vec vector, Vec action);

private:
  static std::string PETScIterativeMethodName(IterativeMethod iterative_method);
};

} // namespace opensn
