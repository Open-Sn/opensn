#pragma once

namespace chi_math
{

template <class VecType, class SolverType>
struct NonLinearSolverContext
{
  virtual ~NonLinearSolverContext() = default;
};

} // namespace chi_math


