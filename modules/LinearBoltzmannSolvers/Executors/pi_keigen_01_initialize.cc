#include "pi_keigen.h"

#include "chi_log_exceptions.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

namespace lbs
{

// ##################################################################
/**Initializer.*/
void
XXPowerIterationKEigen::Initialize()
{
  lbs_solver_.Initialize();

  active_set_source_function_ = lbs_solver_.GetActiveSetSourceFunction();
  primary_ags_solver_ = lbs_solver_.GetPrimaryAGSSolver();

  for (auto& wgs_solver : lbs_solver_.GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<lbs::WGSContext<Mat, Vec, KSP>>(context);

    ChiLogicalErrorIf(not wgs_context, ": Cast failed");

    wgs_context->lhs_src_scope_ =
      wgs_context->lhs_src_scope_ & (~APPLY_WGS_FISSION_SOURCES); // lhs_scope
    wgs_context->rhs_src_scope_ =
      wgs_context->rhs_src_scope_ & (~APPLY_AGS_FISSION_SOURCES); // rhs_scope
  }

  primary_ags_solver_->SetVerbosity(lbs_solver_.Options().verbose_ags_iterations);

  front_wgs_solver_ = lbs_solver_.GetWGSSolvers().at(front_gs_.id_);
  front_wgs_context_ =
    std::dynamic_pointer_cast<lbs::WGSContext<Mat, Vec, KSP>>(front_wgs_solver_->GetContext());

  ChiLogicalErrorIf(not front_wgs_context_, ": Casting failure");

  if (reinit_phi_1_) lbs_solver_.SetPhiVectorScalarValues(phi_old_local_, 1.0);
}

} // namespace lbs