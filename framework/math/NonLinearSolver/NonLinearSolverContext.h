#ifndef CHITECH_CHI_MATH_NL_SOLVER_CONTEXT_H
#define CHITECH_CHI_MATH_NL_SOLVER_CONTEXT_H

namespace chi_math
{

template <class VecType, class SolverType>
struct NonLinearSolverContext
{
  virtual ~NonLinearSolverContext() = default;
};

} // namespace chi_math

#endif // CHITECH_CHI_MATH_NL_SOLVER_CONTEXT_H
