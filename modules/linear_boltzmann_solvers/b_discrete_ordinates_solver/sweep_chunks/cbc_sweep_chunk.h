#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"

namespace opensn
{
class CellMapping;

namespace lbs
{

class CbcSweepChunk : public SweepChunk
{
public:
  CbcSweepChunk(std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const MeshContinuum& grid,
                const SpatialDiscretization& discretization,
                const std::vector<UnitCellMatrices>& unit_cell_matrices,
                std::vector<lbs::CellLBSView>& cell_transport_views,
                const std::vector<double>& source_moments,
                const LBSGroupset& groupset,
                const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                int num_moments,
                int max_num_cell_dofs);

  void SetAngleSet(AngleSet& angle_set) override;

  void SetCell(Cell const* cell_ptr, AngleSet& angle_set) override;

  void Sweep(AngleSet& angle_set) override;

private:
  CBC_FLUDS* fluds_;
  size_t gs_ss_size_;
  size_t gs_ss_begin_;
  int gs_gi_;
  size_t group_stride_;
  size_t group_angle_stride_;
  bool surface_source_active_;

  const Cell* cell_;
  uint64_t cell_local_id_;
  const CellMapping* cell_mapping_;
  CellLBSView* cell_transport_view_;
  size_t cell_num_faces_;
  size_t cell_num_nodes_;

  MatVec3 G_;
  MatDbl M_;
  std::vector<MatDbl> M_surf_;
  std::vector<VecDbl> IntS_shapeI_;
};

} // namespace lbs
} // namespace opensn
