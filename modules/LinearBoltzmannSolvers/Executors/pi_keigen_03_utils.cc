#include "pi_keigen.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs
{

// ##################################################################
/**Combines function calls to set fission source.*/
void XXPowerIterationKEigen::SetLBSFissionSource(const VecDbl& input,
                                                 const bool additive)
{
  if (not additive) chi_math::Set(q_moments_local_, 0.0);
  active_set_source_function_(front_gs_,
                              q_moments_local_,
                              input,
                              APPLY_AGS_FISSION_SOURCES |
                                APPLY_WGS_FISSION_SOURCES);
}

// ##################################################################
/**Combines function calls to set scattering source source.*/
void XXPowerIterationKEigen::SetLBSScatterSource(
  const VecDbl& input, const bool additive, const bool suppress_wg_scat /*=false*/)
{
  if (not additive) chi_math::Set(q_moments_local_, 0.0);
  active_set_source_function_(
    front_gs_,
    q_moments_local_,
    input,
    APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES |
      (suppress_wg_scat ? SUPPRESS_WG_SCATTER : NO_FLAGS_SET));
}

} // namespace lbs