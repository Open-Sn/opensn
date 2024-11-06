// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver_context.h"
#include <string>
#include <utility>
#include <memory>
#include <petscksp.h>

namespace opensn
{

struct LinearSolverContext;

/**
 * Linear solver
 *
 * Wrapper around PETSc KSP
 */
class LinearSolver
{
public:
  enum class IterativeMethod : int
  {
    NONE = 0,
    CLASSIC_RICHARDSON = 1, ///< Classic Richardson (source iteration)
    PETSC_RICHARDSON = 3,   ///< PETSc Richardson iteration
    PETSC_GMRES = 2,        ///< PETSc GMRES iterative algorithm
    PETSC_BICGSTAB = 4,     ///< PETSc BiCGStab iterative algorithm
  };

  struct ToleranceOptions
  {
    double residual_relative = 1.0e-50;
    double residual_absolute = 1.0e-6;
    double residual_divergence = 1.0e6;
    int maximum_iterations = 100;
    int gmres_restart_interval = 100;
    double gmres_breakdown_tolerance = 1.0e6;
  } tolerance_options;

  LinearSolver(IterativeMethod iterative_method, std::shared_ptr<LinearSolverContext> context_ptr);

  virtual ~LinearSolver();

  ToleranceOptions& ToleranceOptions() { return tolerance_options; }
  void ApplyToleranceOptions();

  std::shared_ptr<LinearSolverContext>& GetContext() { return context_ptr_; }

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

protected:
  std::shared_ptr<LinearSolverContext> context_ptr_ = nullptr;
  Mat A_ = nullptr;
  Vec b_ = nullptr;
  Vec x_ = nullptr;
  KSP ksp_ = nullptr;
  int64_t num_local_dofs_ = 0;
  int64_t num_global_dofs_ = 0;

private:
  const IterativeMethod iterative_method_;
  bool use_petsc_ = false;
  bool system_set_ = false;
  bool suppress_kspsolve_ = false;

public:
  static std::string IterativeMethodName(IterativeMethod iterative_method);

private:
  static std::string PETScIterativeMethodName(IterativeMethod iterative_method);
};

} // namespace opensn
