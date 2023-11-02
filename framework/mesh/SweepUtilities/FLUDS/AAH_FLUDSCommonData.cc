#include "AAH_FLUDSCommonData.h"

namespace chi_mesh::sweep_management
{

AAH_FLUDSCommonData::AAH_FLUDSCommonData(
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
  const SPDS& spds,
  const chi_mesh::GridFaceHistogram& grid_face_histogram)
  : FLUDSCommonData(spds, grid_nodal_mappings)
{
  this->InitializeAlphaElements(spds, grid_face_histogram);
  this->InitializeBetaElements(spds);
}

} // namespace chi_mesh::sweep_management