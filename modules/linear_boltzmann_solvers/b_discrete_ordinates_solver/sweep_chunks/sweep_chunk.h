#pragma once

#include "framework/mesh/sweep_utilities/sweep_chunk_base.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_structs.h"

namespace opensn
{
namespace lbs
{

/**Base class for LBS sweepers*/
class SweepChunk : public opensn::SweepChunk
{
public:
  SweepChunk(std::vector<double>& destination_phi,
             std::vector<double>& destination_psi,
             const MeshContinuum& grid,
             const SpatialDiscretization& grid_fe_view,
             const std::vector<UnitCellMatrices>& unit_cell_matrices,
             std::vector<lbs::CellLBSView>& grid_transport_view,
             const std::vector<double>& source_moments,
             const LBSGroupset& groupset,
             const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
             int num_moments,
             int max_num_cell_dofs);

protected:
  const MeshContinuum& grid_;
  const SpatialDiscretization& grid_fe_view_;
  const std::vector<UnitCellMatrices>& unit_cell_matrices_;
  std::vector<lbs::CellLBSView>& grid_transport_view_;
  const std::vector<double>& source_moments_;
  const LBSGroupset& groupset_;
  const std::map<int, std::shared_ptr<MultiGroupXS>>& xs_;
  const int num_moments_;
  const int max_num_cell_dofs_;
  const bool save_angular_flux_;
  const size_t groupset_angle_group_stride_;
  const size_t groupset_group_stride_;
};

} // namespace lbs
} // namespace opensn
