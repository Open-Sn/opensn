#ifndef CHITECH_LBS_NL_KEIGEN_AGS_SOLVER_H
#define CHITECH_LBS_NL_KEIGEN_AGS_SOLVER_H

#include "math/NonLinearSolver/NonLinearSolver.h"
#include "nl_keigen_ags_context.h"

namespace lbs
{

template <class MatType, class VecType, class SolverType>
class NLKEigenvalueAGSSolver
  : public chi_math::NonLinearSolver<MatType, VecType, SolverType>
{
public:
  typedef NLKEigenAGSContext<VecType, SolverType> NLKAGSContext;
  typedef std::shared_ptr<NLKAGSContext> NLKAGSContextPtr;

  explicit NLKEigenvalueAGSSolver(NLKAGSContextPtr nlk_ags_context_ptr)
    : chi_math::NonLinearSolver<MatType, VecType, SolverType>(
        nlk_ags_context_ptr)
  {
  }

  virtual ~NLKEigenvalueAGSSolver() override = default;

protected:
  void PreSetupCallback() override;
  /*void SetOptions();*/
  /*void SetSolverContext();*/
  /*void SetConvergenceTest();*/
  void SetMonitor() override;
  /*void SetPreconditioner();*/

  void SetSystemSize() override;
  void SetSystem() override;
  void SetFunction() override;
  void SetJacobian() override;
  /*void PostSetupCallback();*/
public:
  /*void Setup();*/

protected:
  /*void PreSolveCallback();*/
  void SetInitialGuess() override;
  void PostSolveCallback() override;
  // public:
  //   void Solve();
};

} // namespace lbs

#endif // CHITECH_LBS_NL_KEIGEN_AGS_SOLVER_H
