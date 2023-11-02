#ifndef CHITECH_LBS_ADJOINT_SRC_FUNCTION_H
#define CHITECH_LBS_ADJOINT_SRC_FUNCTION_H

#include "LinearBoltzmannSolvers/A_LBSSolver/SourceFunctions/source_function.h"

namespace lbs
{

/**The adjoint source function removes volumetric fixed source moments
 * as well as point sources, whilst adding volumetric QOI sources.*/
class AdjointSourceFunction : public SourceFunction
{
public:
  explicit AdjointSourceFunction(const LBSSolver& lbs_solver);

  double AddSourceMoments() const override { return 0.0; }

  void AddAdditionalSources(LBSGroupset& groupset,
                            std::vector<double>& destination_q,
                            const std::vector<double>& phi,
                            SourceFlags source_flags) override
  {
    // Inhibit -> AddPointSources
    // Add     -> AddVolumetricQOISources
    AddVolumetricQOISources(groupset, destination_q, phi, source_flags);
  }

  void AddVolumetricQOISources(LBSGroupset& groupset,
                               std::vector<double>& destination_q,
                               const std::vector<double>& phi,
                               SourceFlags source_flags);
};

} // namespace lbs

#endif // CHITECH_LBS_ADJOINT_SRC_FUNCTION_H
