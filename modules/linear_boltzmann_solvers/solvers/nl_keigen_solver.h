// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/nonlinear_keigen_ags_solver.h"
#include <petscsnes.h>

namespace opensn
{

class NonLinearKEigenSolver : public Solver
{
public:
  explicit NonLinearKEigenSolver(const InputParameters& params);

  void Initialize() override;
  void Execute() override;
  /// Return the current k-eigenvalue
  double GetEigenvalue() const;

private:
  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;
  std::shared_ptr<NLKEigenAGSContext> nl_context_;
  NLKEigenvalueAGSSolver nl_solver_;

  bool reset_phi0_;
  int num_initial_power_its_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<NonLinearKEigenSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
