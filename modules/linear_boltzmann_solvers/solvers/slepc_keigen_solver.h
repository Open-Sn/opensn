// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/slepc_linear_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"

namespace opensn
{

class SLEPcKEigenSolver : public Solver
{
private:
  std::shared_ptr<LBSProblem> lbs_problem_;
  std::shared_ptr<SLEPcLinearKEigenContext> context_;
  SLEPcLinearKEigenSolver solver_;

  bool reset_phi0_;

public:
  explicit SLEPcKEigenSolver(const InputParameters& params);

  void Initialize() override;
  void Execute() override;

  /// Return the k-eigenvalue
  double GetEigenvalue() const;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<SLEPcKEigenSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
