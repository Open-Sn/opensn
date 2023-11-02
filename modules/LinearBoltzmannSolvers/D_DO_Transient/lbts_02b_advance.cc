#include "lbts_transient_solver.h"

//###################################################################
/**Advances time values.*/
void
lbs::DiscOrdTransientSolver::Advance()
{
  time_ += dt_;
  phi_prev_local_ = phi_new_local_;
  psi_prev_local_ = psi_new_local_;
  if (options_.use_precursors) precursor_prev_local_ = precursor_new_local_;
}
