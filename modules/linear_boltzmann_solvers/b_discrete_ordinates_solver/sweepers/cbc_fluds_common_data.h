#pragma once

#include "framework/mesh/sweep_utilities/fluds/fluds_common_data.h"

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
