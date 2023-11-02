#include "lbs_solver.h"

//###################################################################
/** Computes the number of moments for the given mesher types*/
void lbs::LBSSolver::ComputeNumberOfMoments()
{
  for (size_t gs = 1; gs < groupsets_.size(); ++gs)
    if (groupsets_[gs].quadrature_->GetMomentToHarmonicsIndexMap()
        != groupsets_[0].quadrature_->GetMomentToHarmonicsIndexMap())
      throw std::logic_error(
        "LinearBoltzmann::SteadyStateSolver::ComputeNumberOfMoments : "
        "Moment-to-Harmonics mapping differs between "
        "groupsets_, which is not allowed.");

  num_moments_ =
    (int)groupsets_.front().quadrature_->GetMomentToHarmonicsIndexMap().size();

  if (num_moments_ == 0)
    throw std::logic_error(
      "LinearBoltzmann::SteadyStateSolver::ComputeNumberOfMoments : "
      "unable to infer number of moments from angular "
      "quadrature.");
}

