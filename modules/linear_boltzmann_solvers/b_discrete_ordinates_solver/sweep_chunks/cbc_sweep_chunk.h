#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"

namespace opensn
{
namespace lbs
{

class CBC_FLUDS;
class CBC_ASynchronousCommunicator;

struct CBC_SweepDependencyInterface : public SweepDependencyInterface
{
  CBC_FLUDS* fluds_ = nullptr;
  const Cell* neighbor_cell_ptr_ = nullptr;

  /**Upwind angular flux*/
  const std::vector<double>* psi_upwnd_data_block_ = nullptr;
  const double* psi_local_face_upwnd_data_ = nullptr;
  /**Downwind angular flux*/
  std::vector<double>* psi_dnwnd_data_ = nullptr;

  size_t group_stride_;
  size_t group_angle_stride_;

  const CellLBSView* cell_transport_view_;

  // Set using SetupIncomingFace
  const FaceNodalMapping* face_nodal_mapping_ = nullptr;

  const double* GetUpwindPsi(int face_node_local_idx) const override;
  double* GetDownwindPsi(int face_node_local_idx) const override;
  void SetupIncomingFace(int face_id,
                         size_t num_face_nodes,
                         uint64_t neighbor_id,
                         bool on_local_face,
                         bool on_boundary) override;
  void SetupOutgoingFace(int face_id,
                         size_t num_face_nodes,
                         uint64_t neighbor_id,
                         bool on_local_face,
                         bool on_boundary,
                         int locality) override;
};

class CBC_SweepChunk : public SweepChunk
{
public:
  CBC_SweepChunk(std::vector<double>& destination_phi,
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

  void SetCells(const std::vector<const Cell*>& cell_ptrs) override;

  void Sweep(AngleSet& angle_set) override;

protected:
  CBC_SweepDependencyInterface sweep_dependency_interface_;
  Cell const* cell_ptr_ = nullptr;
  uint64_t cell_local_id_ = 0;

  std::vector<const Cell*> cell_ptrs_;
};

} // namespace lbs
} // namespace opensn
