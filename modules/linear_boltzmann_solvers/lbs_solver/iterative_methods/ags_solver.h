// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/utils/timer.h"
#include <vector>

namespace opensn
{

/**Solver for Across-Groupset (AGS) solves.*/
class AGSSolver
{
public:
  AGSSolver(LBSSolver& lbs_solver, std::vector<std::shared_ptr<LinearSolver>> wgs_solvers)
    : lbs_solver_(lbs_solver),
      wgs_solvers_(std::move(wgs_solvers)),
      phi_old_(lbs_solver.PhiOldLocal().size()),
      max_iterations_(100),
      tolerance_(1.0e-6),
      verbose_(true)
  {
  }

  ~AGSSolver() {}

  void Solve();

  bool IsVerbose() const { return verbose_; }

  void Verbosity(bool verbose) { verbose_ = verbose; }

  double GetTolerance() { return tolerance_; }

  void Tolerance(double tolerance) { tolerance_ = tolerance; }

  int GetMaxIterations() { return max_iterations_; }

  void MaxIterations(int max_iterations) { max_iterations_ = max_iterations; }

private:
  LBSSolver& lbs_solver_;
  std::vector<std::shared_ptr<LinearSolver>> wgs_solvers_;
  std::vector<double> phi_old_;
  int max_iterations_;
  double tolerance_;
  bool verbose_;

  double ComputeL2PhiChange() const;

  double ComputePointwisePhiChange() const;
};

} // namespace opensn
