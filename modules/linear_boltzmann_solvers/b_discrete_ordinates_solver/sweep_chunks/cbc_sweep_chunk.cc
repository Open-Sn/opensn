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
               max_num_cell_dofs)
{
  sweep_dependency_interface_.groupset_angle_group_stride_ = groupset_angle_group_stride_;
  sweep_dependency_interface_.groupset_group_stride_ = groupset_group_stride_;
}

void
CBC_SweepChunk::SetAngleSet(AngleSet& angle_set)
{
  sweep_dependency_interface_.fluds_ = &dynamic_cast<CBC_FLUDS&>(angle_set.GetFLUDS());

  const SubSetInfo& grp_ss_info = groupset_.grp_subset_infos_[angle_set.GetRefGroupSubset()];

  gs_ss_size_ = grp_ss_info.ss_size;
  gs_ss_begin_ = grp_ss_info.ss_begin;
  gs_gi_ = groupset_.groups_[gs_ss_begin_].id_;

  sweep_dependency_interface_.angle_set_ = &angle_set;
  sweep_dependency_interface_.surface_source_active_ = IsSurfaceSourceActive();
  sweep_dependency_interface_.gs_ss_begin_ = gs_ss_begin_;
  sweep_dependency_interface_.gs_gi_ = gs_gi_;

  sweep_dependency_interface_.group_stride_ = angle_set.GetNumGroups();
  sweep_dependency_interface_.group_angle_stride_ = angle_set.GetNumGroups() * angle_set.GetNumAngles();
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

  sweep_dependency_interface_.cell_transport_view_ = cell_transport_view_;
}

void
CBC_SweepChunk::SetCells(const std::vector<const Cell*>& cell_ptrs)
{
  cell_ptrs_ = cell_ptrs;
}

void
CBC_SweepChunk::Sweep(AngleSet& angle_set)
{
  std::vector<std::vector<double>> Amat(max_num_cell_dofs_, std::vector<double>(max_num_cell_dofs_));
  std::vector<std::vector<double>> Atemp(max_num_cell_dofs_, std::vector<double>(max_num_cell_dofs_));
  std::vector<std::vector<double>> b(groupset_.groups_.size(), std::vector<double>(max_num_cell_dofs_));
  std::vector<double> source(max_num_cell_dofs_);

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
      b[gsg].assign(cell_num_nodes_, 0.0);

    const auto& G = *G_;
    for (int i = 0; i < cell_num_nodes_; ++i)
      for (int j = 0; j < cell_num_nodes_; ++j)
        Amat[i][j] = omega_.Dot(G[i][j]);

    // Update face orientations
    face_mu_values_.assign(cell_num_faces_, 0.0);
    for (int f = 0; f < cell_num_faces_; ++f)
      face_mu_values_[f] = omega_.Dot(cell_->faces_[f].normal_);

    // Surface integrals
    for (int f = 0; f < cell_num_faces_; ++f)
    {
      const auto& face = cell_->faces_[f];

      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const bool local = cell_transport_view_->IsFaceLocal(f);
      const bool boundary = not face.has_neighbor_;

      sweep_dependency_interface_.SetupIncomingFace(
        f, cell_mapping_->NumFaceNodes(f), face.neighbor_id_, local, boundary);

      // IntSf_mu_psi_Mij_dA
      const size_t cf = sweep_dependency_interface_.current_face_idx_;
      const auto& M_surf_f = (*M_surf_)[cf];
      const double mu = face_mu_values_[cf];
      const size_t num_face_nodes = sweep_dependency_interface_.num_face_nodes_;
      for (int fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping_->MapFaceNode(cf, fi);
        for (int fj = 0; fj < num_face_nodes; ++fj)
        {
          const int j = cell_mapping_->MapFaceNode(cf, fj);
          const double* psi = sweep_dependency_interface_.GetUpwindPsi(fj);
          const double mu_Nij = -mu * M_surf_f[i][j];
          Amat[i][j] += mu_Nij;

          if (psi == nullptr)
            continue;

          for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
            b[gsg][i] += psi[gsg] * mu_Nij;
        } // for face node j
      }   // for face node i
    }     // for f

    // Looping over groups, assembling mass terms
    for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
    {
      g_ = gs_gi_ + gsg;
      gsg_ = gsg;
      sigma_tg_ = sigma_t[g_];

      const auto& M = *M_;
      const auto& m2d_op = groupset_.quadrature_->GetMomentToDiscreteOperator();

      // Contribute source moments q = M_n^T * q_moms
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        double temp_src = 0.0;
        for (int m = 0; m < num_moments_; ++m)
        {
          const size_t ir = cell_transport_view_->MapDOF(i, m, static_cast<int>(g_));
          temp_src += m2d_op[m][direction_num_] * q_moments_[ir];
        } // for m
        source[i] = temp_src;
      } // for i

      // Mass Matrix and Source
      // Atemp  = Amat + sigma_tgr * M
      // b     += M * q
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        double temp = 0.0;
        for (int j = 0; j < cell_num_nodes_; ++j)
        {
          const double Mij = M[i][j];
          Atemp[i][j] = Amat[i][j] + Mij * sigma_tg_;
          temp += Mij * source[j];
        } // for j
        b[gsg_][i] += temp;
      } // for i

      // Solve system
      GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes_));
    }

    // Flux updates
    const auto& d2m_op = groupset_.quadrature_->GetDiscreteToMomentOperator();
    auto& output_phi = GetDestinationPhi();
    for (int m = 0; m < num_moments_; ++m)
    {
      const double wn_d2m = d2m_op[m][direction_num_];
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        const size_t ir = cell_transport_view_->MapDOF(i, m, gs_gi_);
        for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
          output_phi[ir + gsg] += wn_d2m * b[gsg][i];
      }
    }

    if (save_angular_flux_)
    {
      auto& output_psi = GetDestinationPsi();
      double* cell_psi_data =
        &output_psi[grid_fe_view_.MapDOFLocal(*cell_, 0, groupset_.psi_uk_man_, 0, 0)];

      for (size_t i = 0; i < cell_num_nodes_; ++i)
      {
        const size_t imap =
          i * groupset_angle_group_stride_ + direction_num_ * groupset_group_stride_ + gs_ss_begin_;
        for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
          cell_psi_data[imap + gsg] = b[gsg][i];
      } // for i
    }

    // Perform outgoing surface operations
    for (int f = 0; f < cell_num_faces_; ++f)
    {
      if (face_orientations[f] != FaceOrientation::OUTGOING)
        continue;

      // Set flags and counters
      const auto& face = cell_->faces_[f];
      const bool local = cell_transport_view_->IsFaceLocal(f);
      const bool boundary = not face.has_neighbor_;
      const int locality = cell_transport_view_->FaceLocality(f);

      sweep_dependency_interface_.SetupOutgoingFace(
        f, cell_mapping_->NumFaceNodes(f), face.neighbor_id_, local, boundary, locality);

      const auto& IntF_shapeI = (*IntS_shapeI_)[f];
      const double mu = face_mu_values_[f];
      const double wt = direction_qweight_;

      const bool on_boundary = sweep_dependency_interface_.on_boundary_;
      const bool is_reflecting_boundary = sweep_dependency_interface_.is_reflecting_bndry_;

      const size_t num_face_nodes = cell_mapping_->NumFaceNodes(f);
      for (int fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping_->MapFaceNode(f, fi);

        double* psi = sweep_dependency_interface_.GetDownwindPsi(fi);

        if (psi != nullptr)
          if (not on_boundary or is_reflecting_boundary)
            for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
              psi[gsg] = b[gsg][i];
        if (on_boundary and not is_reflecting_boundary)
          for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
            cell_transport_view_->AddOutflow(gs_gi_ + gsg, wt * mu * b[gsg][i] * IntF_shapeI[i]);

      } // for fi
    } // for face

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

  if (on_local_face_)
    psi = nullptr; // We don't write local face outputs
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
