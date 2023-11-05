#include "opensn/modules/LinearBoltzmannSolvers/B_DiscreteOrdinatesSolver/Sweepers/CBC_FLUDSCommonData.h"

#include "opensn/framework/mesh/SweepUtilities/SPDS/SPDS.h"
#include "opensn/framework/mesh/MeshContinuum/chi_meshcontinuum.h"

namespace lbs
{

CBC_FLUDSCommonData::CBC_FLUDSCommonData(
  const chi_mesh::sweep_management::SPDS& spds,
  const std::vector<chi_mesh::sweep_management::CellFaceNodalMapping>& grid_nodal_mappings)
  : chi_mesh::sweep_management::FLUDSCommonData(spds, grid_nodal_mappings)
{
}

} // namespace lbs
