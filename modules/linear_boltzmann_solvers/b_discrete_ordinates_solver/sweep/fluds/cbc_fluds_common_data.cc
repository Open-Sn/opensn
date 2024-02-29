#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{
namespace lbs
{

CBC_FLUDSCommonData::CBC_FLUDSCommonData(
  const SPDS& spds, const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
  : FLUDSCommonData(spds, grid_nodal_mappings)
{
}

} // namespace lbs
} // namespace opensn
