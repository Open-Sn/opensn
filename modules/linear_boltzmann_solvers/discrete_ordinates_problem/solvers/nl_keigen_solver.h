// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/nonlinear_keigen_ags_solver.h"
#include <petscsnes.h>

namespace opensn
{

class NonLinearKEigenSolver : public DiscreteOrdinatesSolver
{
public:
  explicit NonLinearKEigenSolver(const InputParameters& params);

  /// Return the current k-eigenvalue
  double GetEigenvalue() const;

  BalanceTable ComputeBalanceTable() const;

protected:
  void CheckRequirements(std::vector<std::string>& errors) const override;
  void InitializeSolver() override;
  void ExecuteSolver() override;

private:
  std::shared_ptr<NLKEigenAGSContext> nl_context_;
  NLKEigenvalueAGSSolver nl_solver_;

  bool reset_phi0_;
  unsigned int num_initial_power_its_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<NonLinearKEigenSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
