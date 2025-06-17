// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "physics/solvers/solver.h"
#include "physics/problems/linear_boltzmann/lbs_problem/lbs_problem.h"
#include "physics/solvers/iterative_methods/nl_keigen_ags_solver.h"
#include <petscsnes.h>

namespace opensn
{

class NonLinearKEigenSolver : public Solver
{
private:
  std::shared_ptr<LBSProblem> lbs_problem_;
  std::shared_ptr<NLKEigenAGSContext> nl_context_;
  NLKEigenvalueAGSSolver nl_solver_;

  bool reset_phi0_;
  int num_initial_power_its_;

public:
  explicit NonLinearKEigenSolver(const InputParameters& params);

  void Initialize() override;
  void Execute() override;
  /// Return the current k-eigenvalue
  double GetEigenvalue() const;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<NonLinearKEigenSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
