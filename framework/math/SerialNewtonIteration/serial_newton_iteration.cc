#include "framework/math/SerialNewtonIteration/serial_newton_iteration.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include <iomanip>

VecDbl
chi_math::NewtonIteration(const NonLinearFunction& non_linear_function,
                          const VecDbl& x_0,
                          const unsigned int max_iters,
                          const double epsilon,
                          const bool verbose /*=false*/)
{
  // Verbose printing lambda
  auto PrintIterationInfo =
    [](unsigned int i, const VecDbl& x_i, const VecDbl& F_x_i, double L2_norm_F_x_i)
  {
    std::stringstream output;
    output << "Iteration " << std::setw(3) << i << ": x_i=";
    for (auto value : x_i)
      output << std::showpos << std::scientific << std::setprecision(3) << value << " ";
    output << "F_x_i=";
    for (auto value : F_x_i)
      output << std::showpos << std::scientific << std::setprecision(3) << value << " ";
    output << "L2_norm_F_x_i=" << L2_norm_F_x_i;

    Chi::log.Log() << output.str();
  };

  // Declare and init variables
  VecDbl x_i = x_0;
  VecDbl F_x_i = non_linear_function.F(x_i);
  MatDbl J_x_i_inv = chi_math::Inverse(non_linear_function.J(x_i));

  double L2_norm_F_x_i = chi_math::Vec2Norm(F_x_i);

  if (verbose) PrintIterationInfo(0, x_i, F_x_i, L2_norm_F_x_i);

  // Perform iterations
  unsigned int i = 0;
  while (L2_norm_F_x_i >= epsilon and i < max_iters)
  {
    ++i;
    x_i = x_i - MatMul(J_x_i_inv, F_x_i);

    F_x_i = non_linear_function.F(x_i);
    J_x_i_inv = chi_math::Inverse(non_linear_function.J(x_i));

    L2_norm_F_x_i = chi_math::Vec2Norm(F_x_i);

    if (verbose) PrintIterationInfo(i, x_i, F_x_i, L2_norm_F_x_i);
  }

  return x_i;
}
