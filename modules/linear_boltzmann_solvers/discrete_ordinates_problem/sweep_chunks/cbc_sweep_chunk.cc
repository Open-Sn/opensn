// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/cell/cell.h"
#include "framework/logging/log.h"
#include "caliper/cali.h"

namespace opensn
{

CBCSweepChunk::CBCSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
  : SweepChunk(problem.GetPhiNewLocal(),
               problem.GetPsiNewLocal()[groupset.id],
               problem.GetGrid(),
               problem.GetSpatialDiscretization(),
               problem.GetUnitCellMatrices(),
               problem.GetCellTransportViews(),
               problem.GetDensitiesLocal(),
               problem.GetQMomentsLocal(),
               groupset,
               problem.GetBlockID2XSMap(),
               problem.GetNumMoments(),
               problem.GetMaxCellDOFCount(),
               problem.GetMinCellDOFCount()),
    fluds_(nullptr),
    gs_size_(0),
    gs_gi_(0),
    num_angles_in_as_(0),
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

  gs_size_ = groupset_.GetNumGroups();
  gs_gi_ = groupset_.first_group;

  surface_source_active_ = IsSurfaceSourceActive();
  num_angles_in_as_ = angle_set.GetNumAngles();
  group_stride_ = angle_set.GetNumGroups();
  group_angle_stride_ = group_stride_ * num_angles_in_as_;
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
  std::vector<Vector<double>> b(gs_size_, Vector<double>(max_num_cell_dofs_));
  std::vector<double> source(max_num_cell_dofs_);

  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[cell_local_id_];
  std::vector<double> face_mu_values(cell_num_faces_);

  const auto& rho = densities_[cell_local_id_];
  const auto& sigma_t = xs_.at(cell_->block_id)->GetSigmaTotal();

  // as = angle set
  // ss = subset
  const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();

  for (size_t as_ss_idx = 0; as_ss_idx < num_angles_in_as_; ++as_ss_idx)
  {
    auto direction_num = as_angle_indices[as_ss_idx];
    auto omega = groupset_.quadrature->omegas[direction_num];
    auto wt = groupset_.quadrature->weights[direction_num];

    // Reset right-hand side
    for (size_t gsg = 0; gsg < gs_size_; ++gsg)
      for (size_t i = 0; i < cell_num_nodes_; ++i)
        b[gsg](i) = 0.0;

    for (size_t i = 0; i < cell_num_nodes_; ++i)
      for (size_t j = 0; j < cell_num_nodes_; ++j)
        Amat(i, j) = omega.Dot(G_(i, j));

    // Update face orientations
    for (size_t f = 0; f < cell_num_faces_; ++f)
      face_mu_values[f] = omega.Dot(cell_->faces[f].normal);

    // Surface integrals
    for (size_t f = 0; f < cell_num_faces_; ++f)
    {
      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const auto& face = cell_->faces[f];
      const bool is_local_face = cell_transport_view_->IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      const auto* face_nodal_mapping =
        &fluds_->GetCommonData().GetFaceNodalMapping(cell_local_id_, f);

      // IntSf_mu_psi_Mij_dA
      const size_t num_face_nodes = cell_mapping_->GetNumFaceNodes(f);
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping_->MapFaceNode(f, fi);

        for (size_t fj = 0; fj < num_face_nodes; ++fj)
        {
          const int j = cell_mapping_->MapFaceNode(f, fj);

          const double mu_Nij = -face_mu_values[f] * M_surf_[f](i, j);
          Amat(i, j) += mu_Nij;

          const double* psi = nullptr;

          if (is_local_face)
            psi = fluds_->UpwindPsi(*cell_transport_view_->FaceNeighbor(f),
                                    face_nodal_mapping->cell_node_mapping_[fj],
                                    as_ss_idx);
          else if (not is_boundary_face)
            psi = fluds_->NLUpwindPsi(
              cell_->global_id, f, face_nodal_mapping->face_node_mapping_[fj], as_ss_idx);
          else
            psi = angle_set.PsiBoundary(face.neighbor_id,
                                        direction_num,
                                        cell_local_id_,
                                        f,
                                        fj,
                                        gs_gi_,
                                        surface_source_active_);

          if (psi != nullptr)
            for (size_t gsg = 0; gsg < gs_size_; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
        } // for face node j
      } // for face node i
    } // for f

    // Looping over groups, assembling mass terms
    for (unsigned int gsg = 0; gsg < gs_size_; ++gsg)
    {
      double sigma_tg = rho * sigma_t[gs_gi_ + gsg];

      // Contribute source moments q = M_n^T * q_moms
      for (size_t i = 0; i < cell_num_nodes_; ++i)
      {
        double temp_src = 0.0;
        for (unsigned int m = 0; m < num_moments_; ++m)
        {
          const auto ir = cell_transport_view_->MapDOF(i, m, gs_gi_ + gsg);
          temp_src += m2d_op[direction_num][m] * source_moments_[ir];
        }
        source[i] = temp_src;
      }

      // Mass matrix and source
      // Atemp = Amat + sigma_tgr * M
      // b += M * q
      for (size_t i = 0; i < cell_num_nodes_; ++i)
      {
        double temp = 0.0;
        for (size_t j = 0; j < cell_num_nodes_; ++j)
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
    for (unsigned int m = 0; m < num_moments_; ++m)
    {
      const double wn_d2m = d2m_op[direction_num][m];
      for (size_t i = 0; i < cell_num_nodes_; ++i)
      {
        const auto ir = cell_transport_view_->MapDOF(i, m, gs_gi_);
        for (size_t gsg = 0; gsg < gs_size_; ++gsg)
          destination_phi_[ir + gsg] += wn_d2m * b[gsg](i);
      }
    }

    // If requested, save angular fluxes during sweep
    if (save_angular_flux_)
    {
      double* cell_psi =
        &destination_psi_[discretization_.MapDOFLocal(*cell_, 0, groupset_.psi_uk_man_, 0, 0)];

      for (size_t i = 0; i < cell_num_nodes_; ++i)
      {
        const size_t addr_offset =
          i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;

        for (size_t gsg = 0; gsg < gs_size_; ++gsg)
          cell_psi[addr_offset + gsg] = b[gsg](i);
      }
    }

    // Perform outgoing surface operations
    for (size_t f = 0; f < cell_num_faces_; ++f)
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
      const auto& face_nodal_mapping =
        fluds_->GetCommonData().GetFaceNodalMapping(cell_local_id_, f);
      std::vector<double>* psi_nonlocal_outgoing = nullptr;

      if (not is_boundary_face and not is_local_face)
      {
        auto& async_comm = *angle_set.GetCommunicator();
        const size_t data_size_for_msg = num_face_nodes * group_angle_stride_;
        psi_nonlocal_outgoing =
          &async_comm.InitGetDownwindMessageData(locality,
                                                 face.neighbor_id,
                                                 face_nodal_mapping.associated_face_,
                                                 angle_set.GetID(),
                                                 data_size_for_msg);
      }

      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping_->MapFaceNode(f, fi);

        // Tally outflow for particle balance
        if (is_boundary_face)
        {
          for (size_t gsg = 0; gsg < gs_size_; ++gsg)
            cell_transport_view_->AddOutflow(
              f, gs_gi_ + gsg, wt * face_mu_values[f] * b[gsg](i) * IntF_shapeI(i));
        }

        double* psi = nullptr;

        if (is_local_face)
          psi = fluds_->OutgoingPsi(*cell_, i, as_ss_idx);
        else if (not is_boundary_face)
          psi = fluds_->NLOutgoingPsi(psi_nonlocal_outgoing, fi, as_ss_idx);
        else if (is_reflecting_boundary_face)
          psi = angle_set.PsiReflected(face.neighbor_id, direction_num, cell_local_id_, f, fi);

        // Write the solved angular flux to the determined location
        if (psi != nullptr)
          for (size_t gsg = 0; gsg < gs_size_; ++gsg)
            psi[gsg] = b[gsg](i);
      } // for fi
    } // for face
  } // for angleset/subset
}

} // namespace opensn
