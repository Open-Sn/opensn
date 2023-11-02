#ifndef CHITECH_LINEAR_SOLVER_CONTEXT_H
#define CHITECH_LINEAR_SOLVER_CONTEXT_H

namespace chi_math
{
  enum class ResidualScaleType
  {
    NONE = 0,
    RHS_NORM = 1,
    RHS_PRECONDITIONED_NORM = 2,
    CUSTOM_SCALE = 3
  };
  template<class MatType, class VecType>
  struct LinearSolverContext
  {
    double rhs_norm = 0.0;
    double rhs_preconditioned_norm = 0.0;
    double custom_residual_scale = 1.0;
    ResidualScaleType residual_scale_type = ResidualScaleType::NONE;

    virtual int MatrixAction(MatType& matrix, VecType& vector, VecType& action)
    {return 0;}

    virtual ~LinearSolverContext() = default;
  };
}//namespace chi_math

#endif //CHITECH_LINEAR_SOLVER_CONTEXT_H
