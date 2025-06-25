// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/cell/cell.h"
#include "framework/logging/log.h"
#include "caliper/cali.h"

namespace opensn
{

CBCSweepChunk::CBCSweepChunk(std::vector<double>& destination_phi,
                             std::vector<double>& destination_psi,
                             const std::shared_ptr<MeshContinuum> grid,
                             const SpatialDiscretization& discretization,
                             const std::vector<UnitCellMatrices>& unit_cell_matrices,
                             std::vector<CellLBSView>& cell_transport_views,
                             const std::vector<double>& densities,
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
               densities,
               source_moments,
               groupset,
               xs,
               num_moments,
               max_num_cell_dofs),
    fluds_(nullptr),
    gs_size_(0),
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
{
}

void
CBCSweepChunk::SetAngleSet(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CbcSweepChunk::SetAngleSet");

  fluds_ = &dynamic_cast<CBC_FLUDS&>(angle_set.GetFLUDS());

  gs_size_ = groupset_.groups.size();
  gs_gi_ = groupset_.groups.front().id;

  surface_source_active_ = IsSurfaceSourceActive();
  group_stride_ = angle_set.GetNumGroups();
  group_angle_stride_ = angle_set.GetNumGroups() * angle_set.GetNumAngles();
}

void
CBCSweepChunk::SetCell(const Cell* cell_ptr, AngleSet& angle_set)
{
  cell_ = cell_ptr;
  cell_local_id_ = cell_ptr->local_id;
  cell_mapping_ = &discretization_.GetCellMapping(*cell_);
  cell_transport_view_ = &cell_transport_views_[cell_->local_id];
  cell_num_faces_ = cell_->faces.size();
  cell_num_nodes_ = cell_mapping_->GetNumNodes();

  // Get cell matrices
  G_ = unit_cell_matrices_[cell_local_id_].intV_shapeI_gradshapeJ;
  M_ = unit_cell_matrices_[cell_local_id_].intV_shapeI_shapeJ;
  M_surf_ = unit_cell_matrices_[cell_local_id_].intS_shapeI_shapeJ;
  IntS_shapeI_ = unit_cell_matrices_[cell_local_id_].intS_shapeI;
}

void
CBCSweepChunk::Sweep(AngleSet& angle_set)
{
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();

  DenseMatrix<double> Amat(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Atemp(max_num_cell_dofs_, max_num_cell_dofs_);
  std::vector<Vector<double>> b(groupset_.groups.size(), Vector<double>(max_num_cell_dofs_));
  std::vector<double> source(max_num_cell_dofs_);

  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[cell_local_id_];
  std::vector<double> face_mu_values(cell_num_faces_);

  const auto& rho = densities_[cell_local_id_];
  const auto& sigma_t = xs_.at(cell_->block_id)->GetSigmaTotal();

  // as = angle set
  // ss = subset
  const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();
  for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
  {
    auto direction_num = as_angle_indices[as_ss_idx];
    auto omega = groupset_.quadrature->omegas[direction_num];
    auto wt = groupset_.quadrature->weights[direction_num];

    // Reset right-hand side
    for (int gsg = 0; gsg < gs_size_; ++gsg)
      for (int i = 0; i < cell_num_nodes_; ++i)
        b[gsg](i) = 0.0;

    for (int i = 0; i < cell_num_nodes_; ++i)
      for (int j = 0; j < cell_num_nodes_; ++j)
        Amat(i, j) = omega.Dot(G_(i, j));

    // Update face orientations
    for (int f = 0; f < cell_num_faces_; ++f)
      face_mu_values[f] = omega.Dot(cell_->faces[f].normal);

    // Surface integrals
    for (int f = 0; f < cell_num_faces_; ++f)
    {
      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const auto& face = cell_->faces[f];
      const bool is_local_face = cell_transport_view_->IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      auto face_nodal_mapping = &fluds_->GetCommonData().GetFaceNodalMapping(cell_local_id_, f);

      const std::vector<double>* psi_upwnd_data_block = nullptr;
      const double* psi_local_face_upwnd_data = nullptr;
      if (is_local_face)
      {
        psi_upwnd_data_block = &fluds_->GetLocalUpwindDataBlock();
        psi_local_face_upwnd_data = fluds_->GetLocalCellUpwindPsi(
          *psi_upwnd_data_block, *cell_transport_view_->FaceNeighbor(f));
      }
      else if (not is_boundary_face)
      {
        psi_upwnd_data_block = &fluds_->GetNonLocalUpwindData(cell_->global_id, f);
      }

      // IntSf_mu_psi_Mij_dA
      const size_t num_face_nodes = cell_mapping_->GetNumFaceNodes(f);
      for (int fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping_->MapFaceNode(f, fi);

        for (int fj = 0; fj < num_face_nodes; ++fj)
        {
          const int j = cell_mapping_->MapFaceNode(f, fj);

          const double mu_Nij = -face_mu_values[f] * M_surf_[f](i, j);
          Amat(i, j) += mu_Nij;

          const double* psi = nullptr;
          if (is_local_face)
          {
            assert(psi_local_face_upwnd_data);
            const unsigned int adj_cell_node = face_nodal_mapping->cell_node_mapping_[fj];
            psi = &psi_local_face_upwnd_data[adj_cell_node * groupset_angle_group_stride_ +
                                             direction_num * groupset_group_stride_];
          }
          else if (not is_boundary_face)
          {
            assert(psi_upwnd_data_block);
            const unsigned int adj_face_node = face_nodal_mapping->face_node_mapping_[fj];
            psi = fluds_->GetNonLocalUpwindPsi(*psi_upwnd_data_block, adj_face_node, as_ss_idx);
          }
          else
            psi = angle_set.PsiBoundary(face.neighbor_id,
                                        direction_num,
                                        cell_local_id_,
                                        f,
                                        fj,
                                        gs_gi_,
                                        surface_source_active_);

          if (not psi)
            continue;

          for (int gsg = 0; gsg < gs_size_; ++gsg)
            b[gsg](i) += psi[gsg] * mu_Nij;
        } // for face node j
      } // for face node i
    } // for f

    // Looping over groups, assembling mass terms
    for (int gsg = 0; gsg < gs_size_; ++gsg)
    {
      double sigma_tg = rho * sigma_t[gs_gi_ + gsg];

      // Contribute source moments q = M_n^T * q_moms
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        double temp_src = 0.0;
        for (int m = 0; m < num_moments_; ++m)
        {
          const size_t ir = cell_transport_view_->MapDOF(i, m, static_cast<int>(gs_gi_ + gsg));
          temp_src += m2d_op[m][direction_num] * source_moments_[ir];
        }
        source[i] = temp_src;
      }

      // Mass matrix and source
      // Atemp = Amat + sigma_tgr * M
      // b += M * q
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        double temp = 0.0;
        for (int j = 0; j < cell_num_nodes_; ++j)
        {
          const double Mij = M_(i, j);
          Atemp(i, j) = Amat(i, j) + Mij * sigma_tg;
          temp += Mij * source[j];
        }
        b[gsg](i) += temp;
      }

      // Solve system
      GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes_));
    } // for gsg

    // Update phi
    for (int m = 0; m < num_moments_; ++m)
    {
      const double wn_d2m = d2m_op[m][direction_num];
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        const size_t ir = cell_transport_view_->MapDOF(i, m, gs_gi_);
        for (int gsg = 0; gsg < gs_size_; ++gsg)
          destination_phi_[ir + gsg] += wn_d2m * b[gsg](i);
      }
    }

    // Save angular flux during sweep
    if (save_angular_flux_)
    {
      double* cell_psi_data =
        &destination_psi_[discretization_.MapDOFLocal(*cell_, 0, groupset_.psi_uk_man_, 0, 0)];

      for (size_t i = 0; i < cell_num_nodes_; ++i)
      {
        const size_t imap =
          i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
        for (int gsg = 0; gsg < gs_size_; ++gsg)
          cell_psi_data[imap + gsg] = b[gsg](i);
      }
    }

    // Perform outgoing surface operations
    for (int f = 0; f < cell_num_faces_; ++f)
    {
      if (face_orientations[f] != FaceOrientation::OUTGOING)
        continue;

      const auto& face = cell_->faces[f];
      const bool is_local_face = cell_transport_view_->IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      const bool is_reflecting_boundary_face =
        (is_boundary_face and angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting());
      const auto& IntF_shapeI = IntS_shapeI_[f];

      const int locality = cell_transport_view_->FaceLocality(f);
      const size_t num_face_nodes = cell_mapping_->GetNumFaceNodes(f);
      auto& face_nodal_mapping = fluds_->GetCommonData().GetFaceNodalMapping(cell_local_id_, f);
      std::vector<double>* psi_dnwnd_data = nullptr;
      if (not is_boundary_face and not is_local_face)
      {
        auto& async_comm = *angle_set.GetCommunicator();
        size_t data_size = num_face_nodes * group_angle_stride_;
        psi_dnwnd_data = &async_comm.InitGetDownwindMessageData(locality,
                                                                face.neighbor_id,
                                                                face_nodal_mapping.associated_face_,
                                                                angle_set.GetID(),
                                                                data_size);
      }

      for (int fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping_->MapFaceNode(f, fi);

        if (is_boundary_face)
        {
          for (int gsg = 0; gsg < gs_size_; ++gsg)
            cell_transport_view_->AddOutflow(
              f, gs_gi_ + gsg, wt * face_mu_values[f] * b[gsg](i) * IntF_shapeI(i));
        }

        double* psi = nullptr;
        if (is_local_face)
          psi = nullptr;
        else if (not is_boundary_face)
        {
          assert(psi_dnwnd_data);
          const size_t addr_offset = fi * group_angle_stride_ + as_ss_idx * group_stride_;
          psi = &(*psi_dnwnd_data)[addr_offset];
        }
        else if (is_reflecting_boundary_face)
          psi = angle_set.PsiReflected(face.neighbor_id, direction_num, cell_local_id_, f, fi);
        if (psi)
        {
          if (not is_boundary_face or is_reflecting_boundary_face)
          {
            for (int gsg = 0; gsg < gs_size_; ++gsg)
              psi[gsg] = b[gsg](i);
          }
        }
      } // for fi
    } // for face
  } // for angleset/subset
}

} // namespace opensn
