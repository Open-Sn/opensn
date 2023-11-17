#pragma once

#include "framework/math/linear_solver/linear_solver.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/wgs_context.h"

#include <memory>
#include <vector>
#include <functional>

namespace lbs
{

/**Linear Solver specialization for Within GroupSet (WGS) solves.*/
class WGSLinearSolver : public chi_math::LinearSolver
{
public:
  typedef std::shared_ptr<WGSContext> WGSContextPtr;

  /**
   * Constructor.
   * \param gs_context_ptr Context Pointer to abstract context.
   */
  explicit WGSLinearSolver(WGSContextPtr gs_context_ptr);
  ~WGSLinearSolver() override;

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

  std::vector<double> saved_q_moments_local_;
};

} // namespace lbs
