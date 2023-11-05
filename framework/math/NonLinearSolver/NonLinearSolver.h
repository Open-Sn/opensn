#pragma once

#include "opensn/framework/math/NonLinearSolver/NonLinearSolverContext.h"
#include "opensn/framework/math/NonLinearSolver/NonLinearSolverOptions.h"

#include <string>

#include <memory>
#include <utility>

namespace chi_math
{

/**Implementation of a general non-linear solver.*/
template <class MatType, class VecType, class SolverType>
class NonLinearSolver
{
public:
  typedef NonLinearSolverContext<VecType, SolverType> NLSolverContext;
  typedef std::shared_ptr<NLSolverContext> NLSolverContextPtr;

  explicit NonLinearSolver(
    NLSolverContextPtr context_ptr,
    const chi::InputParameters& params = NonLinearSolverOptions::GetInputParameters())
    : solver_name_(params.GetParamValue<std::string>("name")),
      context_ptr_(context_ptr),
      options_(params)
  {
  }
  virtual ~NonLinearSolver();

  NonLinearSolverOptions& ToleranceOptions() { return options_; }
  void ApplyToleranceOptions();

  NLSolverContextPtr& GetContext() { return context_ptr_; }

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

  NLSolverContextPtr context_ptr_ = nullptr;

  MatType J_;
  MatType P_;
  VecType r_;
  VecType x_;
  SolverType nl_solver_;

  int64_t num_local_dofs_ = 0;
  int64_t num_globl_dofs_ = 0;

  NonLinearSolverOptions options_;

private:
  bool system_set_ = false;
  bool converged_ = false;
  std::string converged_reason_string_;
};

} // namespace chi_math
