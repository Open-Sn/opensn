#pragma once

#include "linear_solver_context.h"

#include <string>
#include <utility>

#include <memory>

namespace chi_math
{
template <class MatType, class VecType>
struct LinearSolverContext;

/**Implementation of a general linear solver.*/
template <class MatType, class VecType, class SolverType>
class LinearSolver
{
public:
  typedef LinearSolverContext<MatType, VecType> LinSolveContext;
  typedef std::shared_ptr<LinSolveContext> LinSolveContextPtr;

protected:
  const std::string solver_name_;
  const std::string iterative_method_;

  LinSolveContextPtr context_ptr_ = nullptr;

  MatType A_;
  VecType b_;
  VecType x_;
  SolverType solver_;

private:
  bool system_set_ = false;
  bool suppress_kspsolve_ = false;

protected:
  int64_t num_local_dofs_ = 0;
  int64_t num_globl_dofs_ = 0;

  struct ToleranceOptions
  {
    double residual_relative = 1.0e-50;
    double residual_absolute = 1.0e-6;
    double residual_divergence = 1.0e6;
    int maximum_iterations = 100;
    int gmres_restart_interval = 100;
    double gmres_breakdown_tolerance = 1.0e6;
  } tolerance_options_;

protected:
  bool IsSystemSet() const { return system_set_; }

public:
  explicit LinearSolver(const std::string& iterative_method, LinSolveContextPtr context_ptr)
    : solver_name_(iterative_method), iterative_method_(iterative_method), context_ptr_(context_ptr)
  {
  }

  explicit LinearSolver(std::string solver_name,
                        std::string iterative_method,
                        LinSolveContextPtr context_ptr)
    : solver_name_(std::move(solver_name)),
      iterative_method_(std::move(iterative_method)),
      context_ptr_(context_ptr)
  {
  }

  virtual ~LinearSolver();

  ToleranceOptions& ToleranceOptions() { return tolerance_options_; }
  void ApplyToleranceOptions();

  LinSolveContextPtr& GetContext() { return context_ptr_; }

  /**Sets a flag to suppress the KSPSolve() method from being called.*/
  void SetKSPSolveSuppressionFlag(bool flag) { suppress_kspsolve_ = flag; }
  bool GetKSPSolveSuppressionFlag() const { return suppress_kspsolve_; }

protected:
  virtual void PreSetupCallback();
  virtual void SetOptions();
  virtual void SetSolverContext();
  virtual void SetConvergenceTest();
  virtual void SetMonitor();
  virtual void SetPreconditioner();

  virtual void SetSystemSize() = 0;
  virtual void SetSystem() = 0;
  virtual void PostSetupCallback();

public:
  virtual void Setup();

protected:
  virtual void PreSolveCallback();
  virtual void SetInitialGuess() = 0;
  virtual void SetRHS() = 0;
  virtual void PostSolveCallback();

public:
  virtual void Solve();
};

} // namespace chi_math


