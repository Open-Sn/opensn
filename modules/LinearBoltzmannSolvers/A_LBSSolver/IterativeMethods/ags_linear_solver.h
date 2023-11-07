#pragma once

#include "framework/math/LinearSolver/linear_solver.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/ags_context.h"

namespace lbs
{

/**Linear Solver specialization for Within GroupSet (WGS) solves.*/
template <class MatType, class VecType, class SolverType>
class AGSLinearSolver : public chi_math::LinearSolver<MatType, VecType, SolverType>
{
protected:
  int groupspan_first_id_ = 0;
  int groupspan_last_id_ = 0;
  bool verbose_ = false;

public:
  typedef std::shared_ptr<AGSContext<MatType, VecType, SolverType>> AGSContextPtr;

  /**Constructor.
   * \param iterative_method string Across Groupset iterative method.
   * \param ags_context_ptr Pointer Pointer to the context to use.
   * \param groupspan_first_id int First group index.
   * \param groupspan_last_id int Last group index.
   * \param verbose bool Flag to enable verbose output.*/
  AGSLinearSolver(std::string iterative_method,
                  AGSContextPtr ags_context_ptr,
                  int groupspan_first_id,
                  int groupspan_last_id,
                  bool verbose = true)
    : chi_math::LinearSolver<MatType, VecType, SolverType>(std::move(iterative_method),
                                                           ags_context_ptr),
      groupspan_first_id_(groupspan_first_id),
      groupspan_last_id_(groupspan_last_id),
      verbose_(verbose)
  {
  }

  int GroupSpanFirstID() const { return groupspan_first_id_; }
  int GroupSpanLastID() const { return groupspan_last_id_; }
  bool IsVerbose() const { return verbose_; }
  void SetVerbosity(bool verbose_y_n) { verbose_ = verbose_y_n; }

protected:
  /*void PreSetupCallback() override;  */ // Customized via context
  /*virtual void SetOptions();         */
  /*void SetSolverContext() override;  */ // Generic
  /*void SetConvergenceTest() override;*/ // Generic
  /*virtual void SetMonitor();         */

  virtual void SetSystemSize() override; // Customized via context
  virtual void SetSystem() override;     // Generic
  void SetPreconditioner() override;     // Customized via context

  /*void PostSetupCallback() override;*/ // Customized via context

public:
  /*virtual void Setup();*/

protected:
  /*void PreSolveCallback() override;*/  // Customized via context
  void SetRHS() override;                // Generic + with context elements
  void SetInitialGuess() override;       // Generic
  /*void PostSolveCallback() override;*/ // Generic + with context elements
public:
  void Solve() override;

public:
  virtual ~AGSLinearSolver() override;
};

} // namespace lbs
