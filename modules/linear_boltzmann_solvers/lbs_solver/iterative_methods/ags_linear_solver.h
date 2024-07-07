// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_context.h"

namespace opensn
{

/**Linear Solver specialization for Within GroupSet (WGS) solves.*/
class AGSLinearSolver
{
private:
  std::shared_ptr<LinearSolverContext> context_ptr_;
  int maximum_iterations_;
  double tolerance_;
  bool verbose_;

public:
  /**Constructor.
   * \param ags_context_ptr Pointer Pointer to the context to use.
   * \param verbose bool Flag to enable verbose output.*/
  AGSLinearSolver(std::shared_ptr<AGSContext> ags_context_ptr)
    : context_ptr_(ags_context_ptr), maximum_iterations_(1), tolerance_(1.0e-6), verbose_(false)
  {
  }

  ~AGSLinearSolver() {}

  void Solve();

  bool IsVerbose() const { return verbose_; }

  void SetVerbosity(bool verbose) { verbose_ = verbose; }

  double GetTolerance() { return tolerance_; }

  void SetTolerance(double tolerance) { tolerance_ = tolerance; }

  int GetMaximumIterations() { return maximum_iterations_; }

  void SetMaximumIterations(int maximum_iterations) { maximum_iterations_ = maximum_iterations; }
};

} // namespace opensn
