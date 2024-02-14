// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/serial_newton_iteration/serial_newton_iteration.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include <iomanip>

namespace opensn
{

std::vector<double>
NewtonIteration(const NonLinearFunction& non_linear_function,
                const std::vector<double>& x_0,
                const unsigned int max_iters,
                const double epsilon,
                const bool verbose)
{
  // Verbose printing lambda
  auto PrintIterationInfo = [](unsigned int i,
                               const std::vector<double>& x_i,
                               const std::vector<double>& F_x_i,
                               double L2_norm_F_x_i)
  {
    std::stringstream output;
    output << "Iteration " << std::setw(3) << i << ": x_i=";
    for (auto value : x_i)
      output << std::showpos << std::scientific << std::setprecision(3) << value << " ";
    output << "F_x_i=";
    for (auto value : F_x_i)
      output << std::showpos << std::scientific << std::setprecision(3) << value << " ";
    output << "L2_norm_F_x_i=" << L2_norm_F_x_i;

    log.Log() << output.str();
  };

  // Declare and init variables
  std::vector<double> x_i = x_0;
  std::vector<double> F_x_i = non_linear_function.F(x_i);
  MatDbl J_x_i_inv = Inverse(non_linear_function.J(x_i));

  double L2_norm_F_x_i = Vec2Norm(F_x_i);

  if (verbose)
    PrintIterationInfo(0, x_i, F_x_i, L2_norm_F_x_i);

  // Perform iterations
  unsigned int i = 0;
  while (L2_norm_F_x_i >= epsilon and i < max_iters)
  {
    ++i;
    x_i = x_i - MatMul(J_x_i_inv, F_x_i);

    F_x_i = non_linear_function.F(x_i);
    J_x_i_inv = Inverse(non_linear_function.J(x_i));

    L2_norm_F_x_i = Vec2Norm(F_x_i);

    if (verbose)
      PrintIterationInfo(i, x_i, F_x_i, L2_norm_F_x_i);
  }

  return x_i;
}

} // namespace opensn
