#pragma once

#include "framework/math/linear_solver/linear_solver.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_context.h"

#include <memory>
#include <vector>
#include <functional>

namespace lbs
{

/**Linear Solver specialization for Within GroupSet (WGS) solves.*/
template <class MatType, class VecType, class SolverType>
class WGSLinearSolver : public chi_math::LinearSolver<MatType, VecType, SolverType>
{
protected:
  std::vector<double> saved_q_moments_local_;

public:
  typedef std::shared_ptr<WGSContext<MatType, VecType, SolverType>> WGSContextPtr;

  /**Constructor.
   * \param gs_context_ptr Context Pointer to abstract context.*/
  explicit WGSLinearSolver(WGSContextPtr gs_context_ptr)
    : chi_math::LinearSolver<MatType, VecType, SolverType>(
        IterativeMethodPETScName(gs_context_ptr->groupset_.iterative_method_), gs_context_ptr)
  {
    auto& groupset = gs_context_ptr->groupset_;
    auto& solver_tol_options = this->ToleranceOptions();
    solver_tol_options.residual_absolute = groupset.residual_tolerance_;
    solver_tol_options.maximum_iterations = groupset.max_iterations_;
    solver_tol_options.gmres_restart_interval = groupset.gmres_restart_intvl_;
  }

protected:
  /// Customized via context
  void PreSetupCallback() override;
  /// Generic
  void SetConvergenceTest() override;

  /// Customized via context
  void SetSystemSize() override;
  /// Generic
  void SetSystem() override;

  /// Customized via context
  void SetPreconditioner() override;

  /// Customized via context
  void PostSetupCallback() override;

protected:
  /// Customized via context
  void PreSolveCallback() override;
  /// Generic + with context elements
  void SetRHS() override;
  /// Generic
  void SetInitialGuess() override;
  /// Generic + with context elements
  void PostSolveCallback() override;

public:
  ~WGSLinearSolver() override;
};

} // namespace lbs
