#ifndef CHI_MATH_SERIAL_NEWTON_ITERATION_H
#define CHI_MATH_SERIAL_NEWTON_ITERATION_H

#include <vector>
#include <cassert>

#include "math/chi_math.h"


namespace chi_math
{
  /**A simple base class for the evaluation of a non-linear function
   * and its Jacobian-matrix.*/
  class NonLinearFunction
  {
  public:
    /**Function evaluation at vector-x.*/
    virtual VecDbl F(const VecDbl& x) const
    {
      VecDbl result(x.size(), 0.0);
      return result;
    }
    /**Jacobian evaluation at vector-x.*/
    virtual MatDbl J(const VecDbl& x) const
    {
      MatDbl result(x.size(), VecDbl(x.size(), 0.0));
      return result;
    }

    virtual ~NonLinearFunction() = default;
  };

  VecDbl NewtonIteration(const NonLinearFunction& non_linear_function,
                         const VecDbl& x_0,
                         unsigned int max_iters,
                         double epsilon,
                         bool verbose=false);

}//namespace chi_math

#endif //CHI_MATH_SERIAL_NEWTON_ITERATION_H