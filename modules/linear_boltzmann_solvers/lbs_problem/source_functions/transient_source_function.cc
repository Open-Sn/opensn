// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/source_functions/transient_source_function.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"

namespace opensn
{

TransientSourceFunction::TransientSourceFunction(const LBSProblem& lbs_problem)
  : SourceFunction(lbs_problem), lbs_problem_(lbs_problem)
{
}

double
TransientSourceFunction::DelayedFission(const PrecursorList& precursors,
                                        const double& rho,
                                        const std::vector<double>& nu_delayed_sigma_f,
                                        const double* phi) const
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

          value += coeff * eff_dt * precursor.fractional_yield * rho * nu_delayed_sigma_f[gp] *
                   phi[gp] / cell_volume_;
        }

  if (apply_wgs_fission_src_)
    for (size_t gp = gs_i_; gp <= gs_f_; ++gp)
      for (const auto& precursor : precursors)
      {
        const double coeff = precursor.emission_spectrum[g_] * precursor.decay_constant /
                             (1.0 + eff_dt * precursor.decay_constant);

        value += coeff * eff_dt * precursor.fractional_yield * rho * nu_delayed_sigma_f[gp] *
                 phi[gp] / cell_volume_;
      }

  return value;
}

} // namespace opensn
