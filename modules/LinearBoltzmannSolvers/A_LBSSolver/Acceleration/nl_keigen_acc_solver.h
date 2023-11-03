#pragma once

#include "math/NonLinearSolver/NonLinearSolver.h"

#include "nl_keigen_acc_context.h"

#include <petscsnes.h>

namespace lbs::acceleration
{

class NLKEigenDiffSolver : public chi_math::NonLinearSolver<Mat, Vec, SNES>
{
public:
  typedef std::shared_ptr<NLKEigenDiffContext> NLKEigenDiffContextPtr;

  explicit NLKEigenDiffSolver(NLKEigenDiffContextPtr nlk_diff_context_ptr)
    : chi_math::NonLinearSolver<Mat, Vec, SNES>(nlk_diff_context_ptr)
  {
  }

  virtual ~NLKEigenDiffSolver() override = default;

protected:
  /*void PreSetupCallback();*/
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

} // namespace lbs::acceleration


