#pragma once

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/SourceFunctions/source_function.h"

#include "framework/math/math_time_stepping.h"

namespace lbs
{

/**A transient source function needs to adjust the AddDelayedFission
 * routine to properly fit with the current timestepping method and timestep.*/
class TransientSourceFunction : public SourceFunction
{
private:
  double& dt_;
  chi_math::SteppingMethod& method_;

public:
  /**Constructor for the transient source function. The only difference
   * as compared to a steady source function is the treatment of delayed
   * fission.*/
  TransientSourceFunction(const LBSSolver& lbs_solver,
                          double& ref_dt,
                          chi_math::SteppingMethod& method);

  double AddDelayedFission(const PrecursorList& precursors,
                           const std::vector<double>& nu_delayed_sigma_f,
                           const double* phi) const override;
};

} // namespace lbs
