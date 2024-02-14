#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"

namespace opensn
{
namespace lbs
{

class AAH_SweepChunk : public SweepChunk
{
public:
  AAH_SweepChunk(const MeshContinuum& grid,
                 const SpatialDiscretization& discretization,
                 const std::vector<UnitCellMatrices>& unit_cell_matrices,
                 std::vector<lbs::CellLBSView>& cell_transport_views,
                 std::vector<double>& destination_phi,
                 std::vector<double>& destination_psi,
                 const std::vector<double>& source_moments,
                 const LBSGroupset& groupset,
                 const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                 int num_moments,
                 int max_num_cell_dofs);

  void Sweep(AngleSet& angle_set) override;
};

} // namespace lbs
} // namespace opensn
