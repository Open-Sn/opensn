#pragma once

#include "framework/mesh/sweep_utilities/sweep_chunk_base.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_structs.h"

namespace opensn
{
class CellMapping;

namespace lbs
{

struct SweepDependencyInterface
{
  size_t groupset_angle_group_stride_;
  size_t groupset_group_stride_;

  AngleSet* angle_set_ = nullptr;
  bool surface_source_active_ = false;

  size_t gs_ss_begin_ = 0;
  int gs_gi_ = 0;

  const Cell* cell_ptr_ = nullptr;
  uint64_t cell_local_id_ = 0;

  size_t angle_set_index_ = 0;
  size_t angle_num_ = 0;

public: // Set using SetupIncomingFace
  int current_face_idx_ = 0;
  size_t num_face_nodes_ = 0;
  uint64_t neighbor_id_ = 0;
  int face_locality_ = 0;

  bool on_local_face_ = false;
  bool on_boundary_ = false;

public:
  bool is_reflecting_bndry_ = false;

  SweepDependencyInterface() = default;

  virtual const double* GetUpwindPsi(int face_node_local_idx) const = 0;
  virtual double* GetDownwindPsi(int face_node_local_idx) const = 0;

  /**Sets data for the current incoming face.*/
  virtual void SetupIncomingFace(
    int face_id, size_t num_face_nodes, uint64_t neighbor_id, bool on_local_face, bool on_boundary);
  /**Sets data for the current outgoing face.*/
  virtual void SetupOutgoingFace(int face_id,
                                 size_t num_face_nodes,
                                 uint64_t neighbor_id,
                                 bool on_local_face,
                                 bool on_boundary,
                                 int locality);

  virtual ~SweepDependencyInterface() = default;
};

/**Base class for LBS sweepers*/
class SweepChunk : public opensn::SweepChunk
{
public:
  SweepChunk(std::vector<double>& destination_phi,
             std::vector<double>& destination_psi,
             const MeshContinuum& grid,
             const SpatialDiscretization& discretization,
             const std::vector<UnitCellMatrices>& unit_cell_matrices,
             std::vector<lbs::CellLBSView>& cell_transport_views,
             const std::vector<double>& source_moments,
             const LBSGroupset& groupset,
             const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
             int num_moments,
             int max_num_cell_dofs,
             std::unique_ptr<SweepDependencyInterface> sweep_dependency_interface_ptr);

protected:
  const MeshContinuum& grid_;
  const SpatialDiscretization& grid_fe_view_;
  const std::vector<UnitCellMatrices>& unit_cell_matrices_;
  std::vector<lbs::CellLBSView>& grid_transport_view_;
  const std::vector<double>& q_moments_;
  const LBSGroupset& groupset_;
  const std::map<int, std::shared_ptr<MultiGroupXS>>& xs_;
  const int num_moments_;
  const bool save_angular_flux_;

  std::unique_ptr<SweepDependencyInterface> sweep_dependency_interface_ptr_;
  SweepDependencyInterface& sweep_dependency_interface_;

  const size_t groupset_angle_group_stride_;
  const size_t groupset_group_stride_;

  // Runtime params
  size_t gs_ss_size_ = 0;
  size_t gs_ss_begin_ = 0;
  int gs_gi_ = 0;

  std::vector<std::vector<double>> Amat_;
  std::vector<std::vector<double>> Atemp_;
  std::vector<double> source_;
  std::vector<std::vector<double>> b_;

  // Cell items
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
  std::vector<double> face_mu_values_;
  size_t direction_num_ = 0;
  Vector3 omega_;
  double direction_qweight_ = 0.0;
  size_t g_ = 0;
  size_t gsg_ = 0;
  double sigma_tg_ = 0.0;

  // 02 operations
  /**Operations when outgoing fluxes are handled including passing
   * face angular fluxes downstream and computing
   * balance parameters (i.e. outflow)
   * */
  virtual void OutgoingSurfaceOperations();
};

} // namespace lbs
} // namespace opensn
