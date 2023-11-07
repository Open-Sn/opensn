#include "modules/linear_boltzmann_solvers/a_lbs_solver/source_functions/transient_source_function.h"

lbs::TransientSourceFunction::TransientSourceFunction(const LBSSolver& lbs_solver,
                                                      double& ref_dt,
                                                      chi_math::SteppingMethod& method)
  : SourceFunction(lbs_solver), dt_(ref_dt), method_(method)
{
}

double
lbs::TransientSourceFunction::AddDelayedFission(const PrecursorList& precursors,
                                                const std::vector<double>& nu_delayed_sigma_f,
                                                const double* phi) const
{
  const auto& BackwardEuler = chi_math::SteppingMethod::IMPLICIT_EULER;
  const auto& CrankNicolson = chi_math::SteppingMethod::CRANK_NICOLSON;

  double theta;
  if (method_ == BackwardEuler) theta = 1.0;
  else if (method_ == CrankNicolson)
    theta = 0.5;
  else
    theta = 0.7;

  const double eff_dt = theta * dt_;

  double value = 0.0;
  if (apply_ags_fission_src_)
    for (size_t gp = first_grp_; gp <= last_grp_; ++gp)
      if (gp < gs_i_ or gp > gs_f_)
        for (const auto& precursor : precursors)
        {
          const double coeff = precursor.emission_spectrum[g_] * precursor.decay_constant /
                               (1.0 + eff_dt * precursor.decay_constant);

          value += coeff * eff_dt * precursor.fractional_yield * nu_delayed_sigma_f[gp] * phi[gp] /
                   cell_volume_;
        }

  if (apply_wgs_fission_src_)
    for (size_t gp = gs_i_; gp <= gs_f_; ++gp)
      for (const auto& precursor : precursors)
      {
        const double coeff = precursor.emission_spectrum[g_] * precursor.decay_constant /
                             (1.0 + eff_dt * precursor.decay_constant);

        value += coeff * eff_dt * precursor.fractional_yield * nu_delayed_sigma_f[gp] * phi[gp] /
                 cell_volume_;
      }

  return value;
}
