#pragma once

#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/wgs_context.h"

namespace lbs
{
class DiffusionDFEMSolver;
}

namespace lbs
{

struct MIPWGSContext2 : public WGSContext
{
  DiffusionDFEMSolver& lbs_mip_ss_solver_;

  MIPWGSContext2(DiffusionDFEMSolver& lbs_mip_ss_solver,
                 LBSGroupset& groupset,
                 const SetSourceFunction& set_source_function,
                 int lhs_scope,
                 int rhs_scope,
                 bool log_info);

  void PreSetupCallback() override;

  void SetPreconditioner(KSP& solver) override;

  std::pair<int64_t, int64_t> SystemSize() override;

  void ApplyInverseTransportOperator(int scope) override;

  void PostSolveCallback() override;
};
} // namespace lbs
