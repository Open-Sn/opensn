#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"

namespace opensn
{
namespace lbs
{

class CBC_FLUDS;
class CBC_ASynchronousCommunicator;

struct CBC_SweepDependencyInterface
{
  CBC_FLUDS* fluds_ = nullptr;

  const Cell* neighbor_cell_ptr_ = nullptr;
  const Cell* cell_ptr_ = nullptr;
  uint64_t cell_local_id_ = 0;
  const CellLBSView* cell_transport_view_;

  int current_face_idx_ = 0;
  size_t num_face_nodes_ = 0;
  uint64_t neighbor_id_ = 0;
  int face_locality_ = 0;
  bool on_local_face_ = false;
  bool on_boundary_ = false;
  bool is_reflecting_bndry_ = false;
  const FaceNodalMapping* face_nodal_mapping_ = nullptr;

  /**Upwind angular flux*/
  const std::vector<double>* psi_upwnd_data_block_ = nullptr;
  const double* psi_local_face_upwnd_data_ = nullptr;

  /**Downwind angular flux*/
  std::vector<double>* psi_dnwnd_data_ = nullptr;

  size_t groupset_angle_group_stride_;
  size_t groupset_group_stride_;
  size_t group_stride_;
  size_t group_angle_stride_;
  size_t gs_ss_begin_ = 0;
  int gs_gi_ = 0;

  AngleSet* angle_set_ = nullptr;
  size_t angle_set_index_ = 0;
  size_t angle_num_ = 0;

  bool surface_source_active_ = false;
  
  const double* GetUpwindPsi(int face_node_local_idx) const;
  
  double* GetDownwindPsi(int face_node_local_idx) const;
  
  void SetupIncomingFace(int face_id,
                         size_t num_face_nodes,
                         uint64_t neighbor_id,
                         bool on_local_face,
                         bool on_boundary);
  
  void SetupOutgoingFace(int face_id,
                         size_t num_face_nodes,
                         uint64_t neighbor_id,
                         bool on_local_face,
                         bool on_boundary,
                         int locality);
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

private:
  size_t gs_ss_size_ = 0;
  size_t gs_ss_begin_ = 0;
  int gs_gi_ = 0;

  uint64_t cell_local_id_ = 0;
  const Cell* cell_ = nullptr;
  const CellMapping* cell_mapping_ = nullptr;
  CellLBSView* cell_transport_view_ = nullptr;
  size_t cell_num_faces_ = 0;
  size_t cell_num_nodes_ = 0;

  const MatVec3* G_ = nullptr;
  const MatDbl* M_ = nullptr;
  const std::vector<MatDbl>* M_surf_ = nullptr;
  const std::vector<VecDbl>* IntS_shapeI_ = nullptr;

  size_t g_ = 0;
  size_t gsg_ = 0;
 
  CBC_SweepDependencyInterface sweep_dependency_interface_;
  Cell const* cell_ptr_ = nullptr;
  std::vector<const Cell*> cell_ptrs_;
};

} // namespace lbs
} // namespace opensn
