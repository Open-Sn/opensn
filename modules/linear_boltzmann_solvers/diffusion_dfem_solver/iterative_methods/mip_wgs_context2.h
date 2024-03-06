#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"

namespace opensn
{
namespace lbs
{
class DiffusionDFEMSolver;

struct MIPWGSContext2 : public WGSContext
{
  DiffusionDFEMSolver& lbs_mip_ss_solver_;

  MIPWGSContext2(DiffusionDFEMSolver& lbs_mip_ss_solver,
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

} // namespace lbs
} // namespace opensn
