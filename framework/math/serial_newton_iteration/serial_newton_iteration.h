// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/math.h"
#include <vector>
#include <cassert>

namespace opensn
{
/// A simple base class for the evaluation of a non-linear function and its Jacobian-matrix.
class NonLinearFunction
{
public:
  /// Function evaluation at vector-x.
  virtual Vector<double> F(const Vector<double>& x) const
  {
    Vector<double> result(x.Rows(), 0.0);
    return result;
  }

  /// Jacobian evaluation at vector-x.
  virtual DenseMatrix<double> J(const Vector<double>& x) const
  {
    DenseMatrix<double> result(x.Rows(), x.Rows(), 0.0);
    return result;
  }

  virtual ~NonLinearFunction() = default;
};

/// Newton iteration.
Vector<double> NewtonIteration(const NonLinearFunction& non_linear_function,
                               const Vector<double>& x_0,
                               unsigned int max_iters,
                               double epsilon,
                               bool verbose = false);

} // namespace opensn
