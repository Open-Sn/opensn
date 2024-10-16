// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"

namespace opensn
{

class DiffusionDFEMSolver : public LBSSolver
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
};

} // namespace opensn
