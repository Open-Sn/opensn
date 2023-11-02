#ifndef CHITECH_LBS_TRANSIENT_SOURCE_FUNCTION_H
#define CHITECH_LBS_TRANSIENT_SOURCE_FUNCTION_H

#include "source_function.h"

#include "math/chi_math_time_stepping.h"

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
  TransientSourceFunction(const LBSSolver& lbs_solver,
                          double& ref_dt,
                          chi_math::SteppingMethod& method);

  double AddDelayedFission(const PrecursorList& precursors,
                           const std::vector<double>& nu_delayed_sigma_f,
                           const double* phi) const override;
};

}//namespace lbs

#endif //CHITECH_LBS_TRANSIENT_SOURCE_FUNCTION_H
