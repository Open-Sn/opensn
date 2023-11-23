#pragma once

#include "framework/math/nonlinear_solver/nonlinear_solver_context.h"
#include "framework/math/nonlinear_solver/nonlinear_solver_options.h"
#include <string>
#include <memory>
#include <utility>
#include <petscsnes.h>

namespace opensn
{

/**Implementation of a general non-linear solver.*/
class NonLinearSolver
{
public:
  explicit NonLinearSolver(
    std::shared_ptr<NonLinearSolverContext> context_ptr,
    const InputParameters& params = NonLinearSolverOptions::GetInputParameters());
  virtual ~NonLinearSolver();

  NonLinearSolverOptions& ToleranceOptions() { return options_; }
  void ApplyToleranceOptions();

  std::shared_ptr<NonLinearSolverContext>& GetContext() { return context_ptr_; }

  bool IsConverged() const { return converged_; }
  std::string GetConvergedReasonString() const;

  virtual void Setup();
  virtual void Solve();

protected:
  bool IsSystemSet() const { return system_set_; }

  // Called in Setup
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

  // Called in Solver
  virtual void PreSolveCallback();
  virtual void SetInitialGuess() = 0;
  virtual void PostSolveCallback();

  const std::string solver_name_;

  std::shared_ptr<NonLinearSolverContext> context_ptr_ = nullptr;

  Mat J_;
  Mat P_;
  Vec r_;
  Vec x_;
  SNES nl_solver_;

  int64_t num_local_dofs_ = 0;
  int64_t num_globl_dofs_ = 0;

  NonLinearSolverOptions options_;

private:
  bool system_set_ = false;
  bool converged_ = false;
  std::string converged_reason_string_;
};

} // namespace opensn
