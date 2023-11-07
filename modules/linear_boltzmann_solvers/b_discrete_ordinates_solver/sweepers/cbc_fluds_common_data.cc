#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweepers/cbc_fluds_common_data.h"

#include "framework/mesh/sweep_utilities/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace lbs
{

CBC_FLUDSCommonData::CBC_FLUDSCommonData(
  const chi_mesh::sweep_management::SPDS& spds,
  const std::vector<chi_mesh::sweep_management::CellFaceNodalMapping>& grid_nodal_mappings)
  : chi_mesh::sweep_management::FLUDSCommonData(spds, grid_nodal_mappings)
{
}

} // namespace lbs
