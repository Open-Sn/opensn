#pragma once

#include "framework/math/nonlinear_solver/nonlinear_solver.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/nl_keigen_ags_context.h"

namespace lbs
{

template <class MatType, class VecType, class SolverType>
class NLKEigenvalueAGSSolver : public chi_math::NonLinearSolver<MatType, VecType, SolverType>
{
public:
  typedef NLKEigenAGSContext<VecType, SolverType> NLKAGSContext;
  typedef std::shared_ptr<NLKAGSContext> NLKAGSContextPtr;

  explicit NLKEigenvalueAGSSolver(NLKAGSContextPtr nlk_ags_context_ptr)
    : chi_math::NonLinearSolver<MatType, VecType, SolverType>(nlk_ags_context_ptr)
  {
  }

  ~NLKEigenvalueAGSSolver() override = default;

protected:
  void PreSetupCallback() override;
  void SetMonitor() override;

  void SetSystemSize() override;
  void SetSystem() override;
  void SetFunction() override;
  void SetJacobian() override;

protected:
  void SetInitialGuess() override;
  void PostSolveCallback() override;
};

} // namespace lbs
