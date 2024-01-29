#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/cbc_sweep_chunk.h"

#include "framework/mesh/cell/cell.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweepers/cbc_fluds.h"

namespace opensn
{
namespace lbs
{

CBC_SweepChunk::CBC_SweepChunk(std::vector<double>& destination_phi,
                               std::vector<double>& destination_psi,
                               const MeshContinuum& grid,
                               const SpatialDiscretization& discretization,
                               const std::vector<UnitCellMatrices>& unit_cell_matrices,
                               std::vector<lbs::CellLBSView>& cell_transport_views,
                               const std::vector<double>& source_moments,
                               const LBSGroupset& groupset,
                               const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                               int num_moments,
                               int max_num_cell_dofs)
  : SweepChunk(destination_phi,
               destination_psi,
               grid,
               discretization,
               unit_cell_matrices,
               cell_transport_views,
               source_moments,
               groupset,
               xs,
               num_moments,
               max_num_cell_dofs,
               std::make_unique<CBC_SweepDependencyInterface>()),
    cbc_sweep_depinterf_(dynamic_cast<CBC_SweepDependencyInterface&>(sweep_dependency_interface_))
{
  // Register kernels
  RegisterKernel("FEMVolumetricGradTerm",
                 std::bind(&SweepChunk::KernelFEMVolumetricGradientTerm, this));
  RegisterKernel("FEMUpwindSurfaceIntegrals",
                 std::bind(&SweepChunk::KernelFEMUpwindSurfaceIntegrals, this));
  RegisterKernel("FEMSSTDMassTerms", std::bind(&SweepChunk::KernelFEMSTDMassTerms, this));
  RegisterKernel("KernelPhiUpdate", std::bind(&SweepChunk::KernelPhiUpdate, this));
  RegisterKernel("KernelPsiUpdate", std::bind(&SweepChunk::KernelPsiUpdate, this));

  // Setup callbacks
  cell_data_callbacks_ = {};

  direction_data_callbacks_and_kernels_ = {Kernel("FEMVolumetricGradTerm")};

  surface_integral_kernels_ = {Kernel("FEMUpwindSurfaceIntegrals")};

  mass_term_kernels_ = {Kernel("FEMSSTDMassTerms")};

  flux_update_kernels_ = {Kernel("KernelPhiUpdate"), Kernel("KernelPsiUpdate")};

  post_cell_dir_sweep_callbacks_ = {};
}

void
CBC_SweepChunk::SetAngleSet(AngleSet& angle_set)
{
  cbc_sweep_depinterf_.fluds_ = &dynamic_cast<CBC_FLUDS&>(angle_set.GetFLUDS());

  const SubSetInfo& grp_ss_info = groupset_.grp_subset_infos_[angle_set.GetRefGroupSubset()];

  gs_ss_size_ = grp_ss_info.ss_size;
  gs_ss_begin_ = grp_ss_info.ss_begin;
  gs_gi_ = groupset_.groups_[gs_ss_begin_].id_;

  sweep_dependency_interface_.angle_set_ = &angle_set;
  sweep_dependency_interface_.surface_source_active_ = IsSurfaceSourceActive();
  sweep_dependency_interface_.gs_ss_begin_ = gs_ss_begin_;
  sweep_dependency_interface_.gs_gi_ = gs_gi_;

  cbc_sweep_depinterf_.group_stride_ = angle_set.GetNumGroups();
  cbc_sweep_depinterf_.group_angle_stride_ = angle_set.GetNumGroups() * angle_set.GetNumAngles();
}

void
CBC_SweepChunk::SetCell(const Cell* cell_ptr, AngleSet& angle_set)
{
  cell_ptr_ = cell_ptr;
  cell_local_id_ = cell_ptr_->local_id_;

  cell_ = &grid_.local_cells[cell_local_id_];
  sweep_dependency_interface_.cell_ptr_ = cell_;
  sweep_dependency_interface_.cell_local_id_ = cell_local_id_;
  cell_mapping_ = &grid_fe_view_.GetCellMapping(*cell_);
  cell_transport_view_ = &grid_transport_view_[cell_->local_id_];

  cell_num_faces_ = cell_->faces_.size();
  cell_num_nodes_ = cell_mapping_->NumNodes();

  // Get Cell matrices
  const auto& fe_intgrl_values = unit_cell_matrices_[cell_local_id_];
  G_ = &fe_intgrl_values.intV_shapeI_gradshapeJ;
  M_ = &fe_intgrl_values.intV_shapeI_shapeJ;
  M_surf_ = &fe_intgrl_values.intS_shapeI_shapeJ;
  IntS_shapeI_ = &fe_intgrl_values.intS_shapeI;

  for (auto& callback : cell_data_callbacks_)
    callback();

  cbc_sweep_depinterf_.cell_transport_view_ = cell_transport_view_;
}

void
CBC_SweepChunk::SetCells(const std::vector<const Cell*>& cell_ptrs)
{
  cell_ptrs_ = cell_ptrs;
}

void
CBC_SweepChunk::Sweep(AngleSet& angle_set)
{
  using FaceOrientation = FaceOrientation;
  const auto& face_orientations = angle_set.GetSPDS().CellFaceOrientations()[cell_local_id_];
  const auto& sigma_t = xs_.at(cell_->material_id_)->SigmaTotal();

  // as = angle set
  // ss = subset
  const std::vector<size_t>& as_angle_indices = angle_set.GetAngleIndices();
  const size_t as_num_angles = as_angle_indices.size();
  for (size_t as_ss_idx = 0; as_ss_idx < as_num_angles; ++as_ss_idx)
  {
    direction_num_ = as_angle_indices[as_ss_idx];
    omega_ = groupset_.quadrature_->omegas_[direction_num_];
    direction_qweight_ = groupset_.quadrature_->weights_[direction_num_];

    sweep_dependency_interface_.angle_set_index_ = as_ss_idx;
    sweep_dependency_interface_.angle_num_ = direction_num_;

    // Reset right-handside
    for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
      b_[gsg].assign(cell_num_nodes_, 0.0);

    ExecuteKernels(direction_data_callbacks_and_kernels_);

    // Update face orientations
    face_mu_values_.assign(cell_num_faces_, 0.0);
    for (int f = 0; f < cell_num_faces_; ++f)
      face_mu_values_[f] = omega_.Dot(cell_->faces_[f].normal_);

    // Surface integrals
    for (int f = 0; f < cell_num_faces_; ++f)
    {
      const auto& face = cell_->faces_[f];

      if (face_orientations[f] != FaceOrientation::INCOMING) continue;

      const bool local = cell_transport_view_->IsFaceLocal(f);
      const bool boundary = not face.has_neighbor_;

      sweep_dependency_interface_.SetupIncomingFace(
        f, cell_mapping_->NumFaceNodes(f), face.neighbor_id_, local, boundary);

      // IntSf_mu_psi_Mij_dA
      ExecuteKernels(surface_integral_kernels_);
    } // for f

    // Looping over groups, assembling mass terms
    for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
    {
      g_ = gs_gi_ + gsg;
      gsg_ = gsg;
      sigma_tg_ = sigma_t[g_];

      ExecuteKernels(mass_term_kernels_);

      // Solve system
      GaussElimination(Atemp_, b_[gsg], static_cast<int>(cell_num_nodes_));
    }

    // Flux updates
    ExecuteKernels(flux_update_kernels_);

    // Perform outgoing surface operations
    for (int f = 0; f < cell_num_faces_; ++f)
    {
      if (face_orientations[f] != FaceOrientation::OUTGOING) continue;

      // Set flags and counters
      const auto& face = cell_->faces_[f];
      const bool local = cell_transport_view_->IsFaceLocal(f);
      const bool boundary = not face.has_neighbor_;
      const int locality = cell_transport_view_->FaceLocality(f);

      sweep_dependency_interface_.SetupOutgoingFace(
        f, cell_mapping_->NumFaceNodes(f), face.neighbor_id_, local, boundary, locality);

      OutgoingSurfaceOperations();
    } // for face

    ExecuteKernels(post_cell_dir_sweep_callbacks_);
  } // for n
}

void
CBC_SweepDependencyInterface::SetupIncomingFace(
  int face_id, size_t num_face_nodes, uint64_t neighbor_id, bool on_local_face, bool on_boundary)
{
  current_face_idx_ = face_id;
  num_face_nodes_ = num_face_nodes;
  neighbor_id_ = neighbor_id;
  on_local_face_ = on_local_face;
  on_boundary_ = on_boundary;

  face_nodal_mapping_ =
    &fluds_->CommonData().GetFaceNodalMapping(cell_local_id_, current_face_idx_);

  if (on_local_face_)
  {
    neighbor_cell_ptr_ = cell_transport_view_->FaceNeighbor(face_id);
    psi_upwnd_data_block_ = &fluds_->GetLocalUpwindDataBlock();
    psi_local_face_upwnd_data_ =
      fluds_->GetLocalCellUpwindPsi(*psi_upwnd_data_block_, *neighbor_cell_ptr_);
  }
  else if (not on_boundary_)
  {
    psi_upwnd_data_block_ =
      &fluds_->GetNonLocalUpwindData(cell_ptr_->global_id_, current_face_idx_);
  }
}

void
CBC_SweepDependencyInterface::SetupOutgoingFace(int face_id,
                                                size_t num_face_nodes,
                                                uint64_t neighbor_id,
                                                bool on_local_face,
                                                bool on_boundary,
                                                int locality)
{
  current_face_idx_ = face_id;
  num_face_nodes_ = num_face_nodes;
  neighbor_id_ = neighbor_id;
  face_locality_ = locality;
  on_local_face_ = on_local_face;
  on_boundary_ = on_boundary;

  face_nodal_mapping_ =
    &fluds_->CommonData().GetFaceNodalMapping(cell_local_id_, current_face_idx_);

  is_reflecting_bndry_ =
    (on_boundary_ and angle_set_->GetBoundaries()[neighbor_id_]->IsReflecting());

  if (not on_local_face_ and not on_boundary_)
  {
    auto& async_comm = *angle_set_->GetCommunicator();

    size_t data_size = num_face_nodes_ * group_angle_stride_;

    psi_dnwnd_data_ = &async_comm.InitGetDownwindMessageData(face_locality_,
                                                             neighbor_id_,
                                                             face_nodal_mapping_->associated_face_,
                                                             angle_set_->GetID(),
                                                             data_size);
  }
}

const double*
CBC_SweepDependencyInterface::GetUpwindPsi(int face_node_local_idx) const
{
  const double* psi;
  if (on_local_face_)
  {
    const unsigned int adj_cell_node = face_nodal_mapping_->cell_node_mapping_[face_node_local_idx];

    return &psi_local_face_upwnd_data_[adj_cell_node * groupset_angle_group_stride_ +
                                       angle_num_ * groupset_group_stride_ + gs_ss_begin_];
  }
  else if (not on_boundary_)
  {
    const unsigned int adj_face_node = face_nodal_mapping_->face_node_mapping_[face_node_local_idx];

    psi = fluds_->GetNonLocalUpwindPsi(*psi_upwnd_data_block_, adj_face_node, angle_set_index_);
  }
  else
    psi = angle_set_->PsiBndry(neighbor_id_,
                               angle_num_,
                               cell_local_id_,
                               current_face_idx_,
                               face_node_local_idx,
                               gs_gi_,
                               gs_ss_begin_,
                               surface_source_active_);

  return psi;
}

double*
CBC_SweepDependencyInterface::GetDownwindPsi(int face_node_local_idx) const
{
  double* psi = nullptr;

  if (on_local_face_) psi = nullptr; // We don't write local face outputs
  else if (not on_boundary_)
  {
    const size_t addr_offset =
      face_node_local_idx * group_angle_stride_ + angle_set_index_ * group_stride_;

    psi = &(*psi_dnwnd_data_)[addr_offset];
  }
  else if (is_reflecting_bndry_)
    psi = angle_set_->ReflectingPsiOutBoundBndry(neighbor_id_,
                                                 angle_num_,
                                                 cell_local_id_,
                                                 current_face_idx_,
                                                 face_node_local_idx,
                                                 gs_ss_begin_);

  return psi;
}

} // namespace lbs
} // namespace opensn
