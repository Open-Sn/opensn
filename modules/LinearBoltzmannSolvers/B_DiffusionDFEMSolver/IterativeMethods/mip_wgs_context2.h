#ifndef CHITECH_MIP_WGS_CONTEXT_H
#define CHITECH_MIP_WGS_CONTEXT_H

#include "A_LBSSolver/IterativeMethods/wgs_context.h"

namespace lbs
{
class DiffusionDFEMSolver;
}

namespace lbs
{

template <class MatType, class VecType, class SolverType>
struct MIPWGSContext2 : public WGSContext<MatType, VecType, SolverType>
{
  DiffusionDFEMSolver& lbs_mip_ss_solver_;

  MIPWGSContext2(DiffusionDFEMSolver& lbs_mip_ss_solver,
                 LBSGroupset& groupset,
                 const SetSourceFunction& set_source_function,
                 int lhs_scope,
                 int rhs_scope,
                 bool log_info)
    : WGSContext<MatType, VecType, SolverType>(
        lbs_mip_ss_solver, groupset, set_source_function, lhs_scope, rhs_scope, log_info),
      lbs_mip_ss_solver_(lbs_mip_ss_solver)
  {
  }

  void PreSetupCallback() override;

  void SetPreconditioner(SolverType& solver) override;

  std::pair<int64_t, int64_t> SystemSize() override;

  void ApplyInverseTransportOperator(int scope) override;

  void PostSolveCallback() override;
};
} // namespace lbs

#endif // CHITECH_MIP_WGS_CONTEXT_H
