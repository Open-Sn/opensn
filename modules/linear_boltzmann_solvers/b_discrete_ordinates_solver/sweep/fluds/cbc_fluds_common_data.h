#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/fluds/fluds_common_data.h"
#include <cinttypes>

namespace opensn
{
namespace lbs
{

class CBC_FLUDSCommonData : public FLUDSCommonData
{
public:
  CBC_FLUDSCommonData(const SPDS& spds,
                      const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);
};

} // namespace lbs
} // namespace opensn
