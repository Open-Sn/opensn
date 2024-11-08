// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"

namespace opensn
{

class DiffusionDFEMSolver;

struct MIPWGSContext : public WGSContext
{
  MIPWGSContext(DiffusionDFEMSolver& solver,
                LBSGroupset& groupset,
                const SetSourceFunction& set_source_function,
                SourceFlags lhs_scope,
                SourceFlags rhs_scope,
                bool log_info);

  void PreSetupCallback() override;

  void SetPreconditioner(KSP& solver) override;

  std::pair<int64_t, int64_t> SystemSize() override;

  void ApplyInverseTransportOperator(SourceFlags scope) override;

  void PostSolveCallback() override;
};

} // namespace opensn
