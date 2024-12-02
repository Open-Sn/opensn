// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/nl_keigen_ags_solver.h"
#include <petscsnes.h>

namespace opensn
{

class NonLinearKEigen : public opensn::Solver
{
private:
  std::shared_ptr<LBSSolver> lbs_solver_;
  std::shared_ptr<NLKEigenAGSContext> nl_context_;
  NLKEigenvalueAGSSolver nl_solver_;

  bool reset_phi0_;
  int num_initial_power_its_;

public:
  explicit NonLinearKEigen(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<NonLinearKEigen> Create(const ParameterBlock& params);
};

} // namespace opensn
