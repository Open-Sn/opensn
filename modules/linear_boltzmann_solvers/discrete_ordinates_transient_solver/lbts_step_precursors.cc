// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_transient_solver/lbts_transient_solver.h"

#if 0
//###################################################################
/**Performs a timestep of the precursors.*/
void
DiscOrdTransientSolver::StepPrecursors()
{
  const auto& BackwardEuler = opensn::SteppingMethod::IMPLICIT_EULER;
  const auto& CrankNicolson = opensn::SteppingMethod::CRANK_NICOLSON;

  double theta;
  if (method == BackwardEuler) theta = 1.0;
  else if (method == CrankNicolson)
    theta = 0.5;
  else
    theta = 0.7;

  const double eff_dt = theta * dt_;

  // Clear destination vector
  precursor_new_local_.assign(precursor_new_local_.size(), 0.0);

  // Loop over local cells
  // Uses phi_new and precursor_prev_local to compute
  // precursor_new_local(theta-flavor)
  for (auto& cell : grid_ptr_->local_cells)
  {
    const auto& fe_values = unit_cell_matrices_[cell.local_id_];
    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const double cell_volume = transport_view.Volume();

    // Obtain xs
    const auto& xs = matid_to_xs_map_.at(cell.material_id_);
    const auto& precursors = xs->Precursors();
    const auto& nu_delayed_sigma_f = xs->NuDelayedSigmaF();

    // Compute delayed fission rate
    double delayed_fission = 0.0;
    for (int i = 0; i < transport_view.NumNodes(); ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double node_V_fraction = fe_values.Vi_vectors[i] / cell_volume;

      for (int g = 0; g < groups_.size(); ++g)
        delayed_fission += nu_delayed_sigma_f[g] * phi_new_local_[uk_map + g] * node_V_fraction;
    }

    // Loop over precursors
    const auto& max_precursors = max_precursors_per_material_;
    for (unsigned int j = 0; j < xs->NumPrecursors(); ++j)
    {
      const size_t dof_map = cell.local_id_ * max_precursors + j;
      const auto& precursor = precursors[j];
      const double coeff = 1.0 / (1.0 + eff_dt * precursor.decay_constant);

      // contribute last time step precursors
      precursor_new_local_[dof_map] = coeff * precursor_prev_local_[dof_map];

      // contribute delayed fission production
      precursor_new_local_[dof_map] +=
        coeff * eff_dt * precursor.fractional_yield * delayed_fission;
    }
  } // for cell

  // Compute t^{n+1} value
  {
    auto& Cj = precursor_new_local_;
    const auto& Cj_prev = precursor_prev_local_;

    const double inv_theta = 1.0 / theta;
    for (size_t i = 0; i < Cj.size(); ++i)
      Cj[i] = inv_theta * (Cj[i] + (theta - 1.0) * Cj_prev[i]);
  }
}
#endif
