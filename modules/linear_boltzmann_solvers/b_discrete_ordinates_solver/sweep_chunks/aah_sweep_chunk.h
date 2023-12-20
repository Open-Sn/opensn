#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"

#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"

namespace opensn
{
namespace lbs
{

/**Simple utility structure for controlling counters and calls
 * to upstream data.*/
struct AAH_SweepDependencyInterface : public SweepDependencyInterface
{
  AAH_FLUDS* fluds_ = nullptr;

  size_t spls_index = 0;

  int in_face_counter = 0;
  int preloc_face_counter = 0;
  int out_face_counter = 0;
  int deploc_face_counter = 0;

  const double* GetUpwindPsi(int face_node_local_idx) const override;
  double* GetDownwindPsi(int face_node_local_idx) const override;
};

/**The new sweep chunk class.*/
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
                 const std::map<int, XSPtr>& xs,
                 int num_moments,
                 int max_num_cell_dofs);

  // 01
  void Sweep(AngleSet& angle_set) override;
};

} // namespace lbs
} // namespace opensn
