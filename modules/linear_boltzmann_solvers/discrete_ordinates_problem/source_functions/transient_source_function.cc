// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/source_functions/transient_source_function.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"

namespace opensn
{

TransientSourceFunction::TransientSourceFunction(const LBSProblem& lbs_problem)
  : SourceFunction(lbs_problem)
{
}

double
TransientSourceFunction::DelayedFission(const PrecursorList& precursors,
                                        const std::vector<double>& nu_delayed_sigma_f,
                                        const double* phi,
                                        std::uint64_t cell_local_id) const
{
  const double eff_dt = lbs_problem_.GetTheta() * lbs_problem_.GetTimeStep();

  double value = 0.0;
  if (apply_ags_fission_src_)
    for (size_t gp = first_grp_; gp <= last_grp_; ++gp)
      if (gp < gs_i_ or gp > gs_f_)
        for (const auto& precursor : precursors)
        {
          const double coeff = precursor.emission_spectrum[g_] * precursor.decay_constant /
                               (1.0 + eff_dt * precursor.decay_constant);

          value += coeff * eff_dt * precursor.fractional_yield * nu_delayed_sigma_f[gp] * phi[gp];
        }

  if (apply_wgs_fission_src_)
    for (size_t gp = gs_i_; gp <= gs_f_; ++gp)
      for (const auto& precursor : precursors)
      {
        const double coeff = precursor.emission_spectrum[g_] * precursor.decay_constant /
                             (1.0 + eff_dt * precursor.decay_constant);

        value += coeff * eff_dt * precursor.fractional_yield * nu_delayed_sigma_f[gp] * phi[gp];
      }

  // Contribution from the decay of precursor inventory accumulated prior to this time step
  // (C_j(t_n)). This term is independent of the flux being solved for, so it belongs with the
  // other fixed (RHS-only) sources: folding it in under the AGS/WGS fission flags above would
  // incorrectly re-add this flux-independent constant on every inner (LHS-operator) iteration.
  if (apply_fixed_src_ and not precursors.empty())
  {
    const auto& precursor_old_local = lbs_problem_.GetPrecursorsOldLocal();
    const auto max_precursors = lbs_problem_.GetMaxPrecursorsPerMaterial();
    const auto cell_base = cell_local_id * max_precursors;
    for (std::size_t j = 0; j < precursors.size(); ++j)
    {
      const auto& precursor = precursors[j];
      const double coeff = precursor.emission_spectrum[g_] * precursor.decay_constant /
                           (1.0 + eff_dt * precursor.decay_constant);
      value += coeff * precursor_old_local[cell_base + j];
    }
  }

  return value;
}

} // namespace opensn
