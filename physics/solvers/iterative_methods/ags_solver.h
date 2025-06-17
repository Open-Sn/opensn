// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "physics/problems/linear_boltzmann/lbs_problem/lbs_problem.h"
#include "framework/utils/timer.h"
#include <vector>

namespace opensn
{

/**Solver for Across-Groupset (AGS) solves.*/
class AGSSolver
{
public:
  AGSSolver(LBSProblem& lbs_problem, std::vector<std::shared_ptr<LinearSolver>> wgs_solvers)
    : lbs_problem_(lbs_problem),
      wgs_solvers_(std::move(wgs_solvers)),
      phi_old_(lbs_problem.GetPhiOldLocal().size()),
      max_iterations_(100),
      tolerance_(1.0e-6),
      verbose_(true)
  {
  }

  ~AGSSolver() {}

  void Solve();

  bool IsVerbose() const { return verbose_; }

  void SetVerbosity(bool verbose) { verbose_ = verbose; }

  double GetTolerance() { return tolerance_; }

  void SetTolerance(double tolerance) { tolerance_ = tolerance; }

  int GetMaxIterations() { return max_iterations_; }

  void SetMaxIterations(int max_iterations) { max_iterations_ = max_iterations; }

private:
  LBSProblem& lbs_problem_;
  std::vector<std::shared_ptr<LinearSolver>> wgs_solvers_;
  std::vector<double> phi_old_;
  int max_iterations_;
  double tolerance_;
  bool verbose_;
};

} // namespace opensn
