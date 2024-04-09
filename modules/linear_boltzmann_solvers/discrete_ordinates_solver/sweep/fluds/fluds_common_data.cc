// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/fluds/fluds_common_data.h"

namespace opensn
{
namespace lbs
{

FLUDSCommonData::FLUDSCommonData(const SPDS& spds,
                                 const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
  : spds_(spds), grid_nodal_mappings_(grid_nodal_mappings)
{
}

const SPDS&
FLUDSCommonData::GetSPDS() const
{
  return spds_;
}

const FaceNodalMapping&
FLUDSCommonData::GetFaceNodalMapping(uint64_t cell_local_id, unsigned int face_id) const
{
  return grid_nodal_mappings_[cell_local_id][face_id];
}

} // namespace lbs
} // namespace opensn
