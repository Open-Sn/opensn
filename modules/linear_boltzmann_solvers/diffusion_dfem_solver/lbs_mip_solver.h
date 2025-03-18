// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"

namespace opensn
{

class DiffusionDFEMSolver : public LBSProblem
{
public:
  std::vector<std::shared_ptr<DiffusionMIPSolver>> gs_mip_solvers;

public:
  explicit DiffusionDFEMSolver(const InputParameters& params);
  ~DiffusionDFEMSolver() override;
  void Initialize() override;
  void InitializeWGSSolvers() override;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<DiffusionDFEMSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
