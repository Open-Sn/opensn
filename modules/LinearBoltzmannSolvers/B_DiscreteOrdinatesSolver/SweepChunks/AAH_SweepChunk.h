#ifndef CHITECH_AAH_SWEEPCHUNK_H
#define CHITECH_AAH_SWEEPCHUNK_H

#include "SweepChunk.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

namespace lbs
{

// ##################################################################
/**Simple utility structure for controlling counters and calls
 * to upstream data.*/
struct AAH_SweepDependencyInterface : public SweepDependencyInterface
{
  chi_mesh::sweep_management::AAH_FLUDS* fluds_ = nullptr;

  size_t spls_index = 0;

  int in_face_counter = 0;
  int preloc_face_counter = 0;
  int out_face_counter = 0;
  int deploc_face_counter = 0;

  const double* GetUpwindPsi(int face_node_local_idx) const override;
  double* GetDownwindPsi(int face_node_local_idx) const override;
};

// ##################################################################
/**The new sweep chunk class.*/
class AAH_SweepChunk : public SweepChunk
{
public:
  AAH_SweepChunk(const chi_mesh::MeshContinuum& grid,
                 const chi_math::SpatialDiscretization& discretization,
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
  void Sweep(chi_mesh::sweep_management::AngleSet& angle_set) override;
};

} // namespace lbs

#endif // CHITECH_AAH_SWEEPCHUNK_H
