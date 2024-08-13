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
  virtual DenseVector<double> F(const DenseVector<double>& x) const
  {
    DenseVector<double> result(x.Rows(), 0.0);
    return result;
  }

  /// Jacobian evaluation at vector-x.
  virtual DenseMatrix<double> J(const DenseVector<double>& x) const
  {
    DenseMatrix<double> result(x.Rows(), x.Rows(), 0.0);
    return result;
  }

  virtual ~NonLinearFunction() = default;
};

/// Newton iteration.
DenseVector<double> NewtonIteration(const NonLinearFunction& non_linear_function,
                                    const DenseVector<double>& x_0,
                                    unsigned int max_iters,
                                    double epsilon,
                                    bool verbose = false);

} // namespace opensn
