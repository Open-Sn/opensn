#include "framework/mesh/SweepUtilities/FLUDS/FLUDSCommonData.h"

namespace chi_mesh::sweep_management
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

} // namespace chi_mesh::sweep_management
