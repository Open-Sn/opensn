#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/cbc_sweep_chunk.h"
#include "framework/mesh/cell/cell.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"

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
               max_num_cell_dofs),
    fluds_(nullptr),
    gs_ss_size_(0),
    gs_ss_begin_(0),
    gs_gi_(0),
    group_stride_(0),
    group_angle_stride_(0),
    surface_source_active_(false),
    cell_(nullptr),
    cell_local_id_(0),
    cell_mapping_(nullptr),
    cell_transport_view_(nullptr),
    cell_num_faces_(0),
    cell_num_nodes_(0)
{}

void
CBC_SweepChunk::SetAngleSet(AngleSet& angle_set)
{
  fluds_ = &dynamic_cast<CBC_FLUDS&>(angle_set.GetFLUDS());

  const SubSetInfo& grp_ss_info = groupset_.grp_subset_infos_[angle_set.GetRefGroupSubset()];
  
  gs_ss_size_ = grp_ss_info.ss_size;
  gs_ss_begin_ = grp_ss_info.ss_begin;
  gs_gi_ = groupset_.groups_[gs_ss_begin_].id_;

  surface_source_active_ = IsSurfaceSourceActive();
  group_stride_ = angle_set.GetNumGroups();
  group_angle_stride_ = angle_set.GetNumGroups() * angle_set.GetNumAngles();
}

void
CBC_SweepChunk::SetCell(const Cell* cell_ptr, AngleSet& angle_set)
{
  cell_ = cell_ptr;
  cell_local_id_ = cell_ptr->local_id_;
  cell_mapping_ = &grid_fe_view_.GetCellMapping(*cell_);
  cell_transport_view_ = &grid_transport_view_[cell_->local_id_];
  cell_num_faces_ = cell_->faces_.size();
  cell_num_nodes_ = cell_mapping_->NumNodes();

  // Get cell matrices
  G_ = unit_cell_matrices_[cell_local_id_].intV_shapeI_gradshapeJ;
  M_ = unit_cell_matrices_[cell_local_id_].intV_shapeI_shapeJ;
  M_surf_ = unit_cell_matrices_[cell_local_id_].intS_shapeI_shapeJ;
  IntS_shapeI_ = unit_cell_matrices_[cell_local_id_].intS_shapeI;
}

void
CBC_SweepChunk::SetCells(const std::vector<const Cell*>& cell_ptrs)
{}

void
CBC_SweepChunk::Sweep(AngleSet& angle_set)
{
  const auto& m2d_op = groupset_.quadrature_->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature_->GetDiscreteToMomentOperator();

  std::vector<std::vector<double>> Amat(max_num_cell_dofs_, std::vector<double>(max_num_cell_dofs_));
  std::vector<std::vector<double>> Atemp(max_num_cell_dofs_, std::vector<double>(max_num_cell_dofs_));
  std::vector<std::vector<double>> b(groupset_.groups_.size(), std::vector<double>(max_num_cell_dofs_));
  std::vector<double> source(max_num_cell_dofs_);

  const auto& face_orientations = angle_set.GetSPDS().CellFaceOrientations()[cell_local_id_];
   std::vector<double> face_mu_values(cell_num_faces_);

  const auto& sigma_t = xs_.at(cell_->material_id_)->SigmaTotal();

  // as = angle set
  // ss = subset
  const std::vector<size_t>& as_angle_indices = angle_set.GetAngleIndices();
  for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
  {
    auto direction_num = as_angle_indices[as_ss_idx];
    auto omega = groupset_.quadrature_->omegas_[direction_num];
    auto wt = groupset_.quadrature_->weights_[direction_num];

    // Reset right-hand side
    for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
      b[gsg].assign(cell_num_nodes_, 0.0);

    for (int i = 0; i < cell_num_nodes_; ++i)
      for (int j = 0; j < cell_num_nodes_; ++j)
        Amat[i][j] = omega.Dot(G_[i][j]);

    // Update face orientations
    for (int f = 0; f < cell_num_faces_; ++f)
      face_mu_values[f] = omega.Dot(cell_->faces_[f].normal_);

    // Surface integrals
    for (int f = 0; f < cell_num_faces_; ++f)
    {
      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const auto& face = cell_->faces_[f];
      const bool local = cell_transport_view_->IsFaceLocal(f);
      const bool boundary = not face.has_neighbor_;
      auto face_nodal_mapping = &fluds_->CommonData().GetFaceNodalMapping(cell_local_id_, f);

      const std::vector<double>* psi_upwnd_data_block = nullptr;
      const double* psi_local_face_upwnd_data = nullptr;
      if (local)
      {
        psi_upwnd_data_block = &fluds_->GetLocalUpwindDataBlock();
        psi_local_face_upwnd_data = fluds_->GetLocalCellUpwindPsi(*psi_upwnd_data_block, *cell_transport_view_->FaceNeighbor(f));
      }
      else if (not boundary)
      {
        psi_upwnd_data_block = &fluds_->GetNonLocalUpwindData(cell_->global_id_, f);
      }

      // IntSf_mu_psi_Mij_dA
      const auto& M_surf_f = M_surf_[f];
      const double mu = face_mu_values[f];
      const size_t num_face_nodes = cell_mapping_->NumFaceNodes(f);
      for (int fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping_->MapFaceNode(f, fi);
        for (int fj = 0; fj < num_face_nodes; ++fj)
        {
          const int j = cell_mapping_->MapFaceNode(f, fj);

          const double* psi = nullptr;
          if (local)
          {
            const unsigned int adj_cell_node = face_nodal_mapping->cell_node_mapping_[fj];
            psi = &psi_local_face_upwnd_data[adj_cell_node*groupset_angle_group_stride_ +
              direction_num*groupset_group_stride_ + gs_ss_begin_];
          }
          else if (not boundary)
          {
            const unsigned int adj_face_node = face_nodal_mapping->face_node_mapping_[fj];
            psi = fluds_->GetNonLocalUpwindPsi(*psi_upwnd_data_block, adj_face_node, as_ss_idx);
          }
          else
            psi = angle_set.PsiBndry(face.neighbor_id_,
                                     direction_num,
                                     cell_local_id_,
                                     f,
                                     fj,
                                     gs_gi_,
                                     gs_ss_begin_,
                                     surface_source_active_);

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
      double sigma_tg = sigma_t[gs_gi_ + gsg];

      // Contribute source moments q = M_n^T * q_moms
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        double temp_src = 0.0;
        for (int m = 0; m < num_moments_; ++m)
        {
          const size_t ir = cell_transport_view_->MapDOF(i, m, static_cast<int>(gs_gi_ + gsg));
          temp_src += m2d_op[m][direction_num] * source_moments_[ir];
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
          const double Mij = M_[i][j];
          Atemp[i][j] = Amat[i][j] + Mij * sigma_tg;
          temp += Mij * source[j];
        } // for j
        b[gsg][i] += temp;
      } // for i

      // Solve system
      GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes_));
    }

    // Flux updates
    auto& output_phi = GetDestinationPhi();
    for (int m = 0; m < num_moments_; ++m)
    {
      const double wn_d2m = d2m_op[m][direction_num];
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
          i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_ + gs_ss_begin_;
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
      const size_t num_face_nodes = cell_mapping_->NumFaceNodes(f);
      auto face_nodal_mapping = &fluds_->CommonData().GetFaceNodalMapping(cell_local_id_, f);
      const bool is_reflecting_boundary =
        (boundary and angle_set.GetBoundaries()[face.neighbor_id_]->IsReflecting());

      std::vector<double>* psi_dnwnd_data = nullptr;
      if (not local and not boundary)
      {
        auto& async_comm = *angle_set.GetCommunicator();
        size_t data_size = num_face_nodes * group_angle_stride_;
        psi_dnwnd_data = &async_comm.InitGetDownwindMessageData(locality,
                                                                face.neighbor_id_,
                                                                face_nodal_mapping->associated_face_,
                                                                angle_set.GetID(),
                                                                data_size);
      }

      const double mu = face_mu_values[f];
      for (int fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping_->MapFaceNode(f, fi);

        double* psi = nullptr;
        if (local)
          psi = nullptr; // We don't write local face outputs
        else if (not boundary)
        {
          const size_t addr_offset = fi * group_angle_stride_ + as_ss_idx * group_stride_;
          psi = &(*psi_dnwnd_data)[addr_offset];
        }
        else if (is_reflecting_boundary)
          psi = angle_set.ReflectingPsiOutBoundBndry(face.neighbor_id_,
                                                     direction_num,
                                                     cell_local_id_,
                                                     f,
                                                     fi,
                                                     gs_ss_begin_);
        if (psi != nullptr)
          if (not boundary or is_reflecting_boundary)
            for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
              psi[gsg] = b[gsg][i];
        if (boundary and not is_reflecting_boundary)
          for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
            cell_transport_view_->AddOutflow(gs_gi_ + gsg, wt * mu * b[gsg][i] * IntS_shapeI_[f][i]);
      } // for fi
    } // for face
  } // for n
}

} // namespace lbs
} // namespace opensn
