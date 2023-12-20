#pragma once

#include "modules/linear_boltzmann_solvers/a_lbs_solver/source_functions/source_function.h"

namespace opensn
{
namespace lbs
{

/**The adjoint source function removes volumetric fixed source moments
 * as well as point sources, whilst adding volumetric QOI sources.*/
class AdjointSourceFunction : public SourceFunction
{
public:
  /**Constructor for an adjoint source function.*/
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

  /**Adds Quantities of Interest to the nodal sources.*/
  void AddVolumetricQOISources(LBSGroupset& groupset,
                               std::vector<double>& destination_q,
                               const std::vector<double>& phi,
                               SourceFlags source_flags);
};

} // namespace lbs
} // namespace opensn
