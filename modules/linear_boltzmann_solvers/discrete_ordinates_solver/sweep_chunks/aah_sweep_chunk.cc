// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/fluds/aah_fluds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"

namespace opensn
{

AahSweepChunk::AahSweepChunk(const MeshContinuum& grid,
                             const SpatialDiscretization& discretization,
                             const std::vector<UnitCellMatrices>& unit_cell_matrices,
                             std::vector<CellLBSView>& cell_transport_views,
                             const std::vector<double>& densities,
                             std::vector<double>& destination_phi,
                             std::vector<double>& destination_psi,
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
               max_num_cell_dofs)
{
}

void
AahSweepChunk::Sweep(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("AahSweepChunk::Sweep");

  const SubSetInfo& grp_ss_info = groupset_.grp_subset_infos[angle_set.GetGroupSubset()];

  auto gs_ss_size = grp_ss_info.ss_size;
  auto gs_ss_begin = grp_ss_info.ss_begin;
  auto gs_gi = groupset_.groups[gs_ss_begin].id;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();

  DenseMatrix<double> Amat(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Atemp(max_num_cell_dofs_, max_num_cell_dofs_);
  std::vector<Vector<double>> b(groupset_.groups.size(), Vector<double>(max_num_cell_dofs_, 0.));
  std::vector<double> source(max_num_cell_dofs_);

  // Loop over each cell
  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.LocalSubgrid();
  const size_t num_spls = spls.size();
  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    auto cell_local_id = spls[spls_index];
    auto& cell = grid_.local_cells[cell_local_id];
    auto& cell_mapping = discretization_.GetCellMapping(cell);
    auto& cell_transport_view = cell_transport_views_[cell_local_id];
    auto cell_num_faces = cell.faces.size();
    auto cell_num_nodes = cell_mapping.NumNodes();

    const auto& face_orientations = spds.CellFaceOrientations()[cell_local_id];
    std::vector<double> face_mu_values(cell_num_faces);

    const auto& rho = densities_[cell.local_id];
    const auto& sigma_t = xs_.at(cell.material_id)->SigmaTotal();

    // Get cell matrices
    const auto& G = unit_cell_matrices_[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = unit_cell_matrices_[cell_local_id].intS_shapeI_shapeJ;

    // Loop over angles in set (as = angleset, ss = subset)
    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;
    const std::vector<size_t>& as_angle_indices = angle_set.GetAngleIndices();
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      auto direction_num = as_angle_indices[as_ss_idx];
      auto omega = groupset_.quadrature->omegas[direction_num];
      auto wt = groupset_.quadrature->weights[direction_num];

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      // Reset right-hand side
      for (int gsg = 0; gsg < gs_ss_size; ++gsg)
        for (int i = 0; i < cell_num_nodes; ++i)
          b[gsg](i) = 0.0;

      for (int i = 0; i < cell_num_nodes; ++i)
        for (int j = 0; j < cell_num_nodes; ++j)
          Amat(i, j) = omega.Dot(G(i, j));

      // Update face orientations
      for (int f = 0; f < cell_num_faces; ++f)
        face_mu_values[f] = omega.Dot(cell.faces[f].normal);

      // Surface integrals
      int in_face_counter = -1;
      for (int f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::INCOMING)
          continue;

        auto& cell_face = cell.faces[f];
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary_face = not cell_face.has_neighbor;

        if (is_local_face)
          ++in_face_counter;
        else if (not is_boundary_face)
          ++preloc_face_counter;

        // IntSf_mu_psi_Mij_dA
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (int fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);

          for (int fj = 0; fj < num_face_nodes; ++fj)
          {
            const int j = cell_mapping.MapFaceNode(f, fj);

            const double mu_Nij = -face_mu_values[f] * M_surf[f](i, j);
            Amat(i, j) += mu_Nij;

            const double* psi;
            if (is_local_face)
              psi = fluds.UpwindPsi(spls_index, in_face_counter, fj, 0, as_ss_idx);
            else if (not is_boundary_face)
              psi = fluds.NLUpwindPsi(preloc_face_counter, fj, 0, as_ss_idx);
            else
              psi = angle_set.PsiBoundary(cell_face.neighbor_id,
                                          direction_num,
                                          cell_local_id,
                                          f,
                                          fj,
                                          gs_gi,
                                          gs_ss_begin,
                                          IsSurfaceSourceActive());

            if (not psi)
              continue;

            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
          } // for face node j
        }   // for face node i
      }     // for f

      // Looping over groups, assembling mass terms
      for (int gsg = 0; gsg < gs_ss_size; ++gsg)
      {
        double sigma_tg = rho * sigma_t[gs_gi + gsg];

        // Contribute source moments q = M_n^T * q_moms
        for (int i = 0; i < cell_num_nodes; ++i)
        {
          double temp_src = 0.0;
          for (int m = 0; m < num_moments_; ++m)
          {
            const size_t ir = cell_transport_view.MapDOF(i, m, static_cast<int>(gs_gi + gsg));
            temp_src += m2d_op[m][direction_num] * source_moments_[ir];
          }
          source[i] = temp_src;
        }

        // Mass matrix and source
        // Atemp = Amat + sigma_tgr * M
        // b += M * q
        for (int i = 0; i < cell_num_nodes; ++i)
        {
          double temp = 0.0;
          for (int j = 0; j < cell_num_nodes; ++j)
          {
            auto Mij = M(i, j);
            Atemp(i, j) = Amat(i, j) + Mij * sigma_tg;
            temp += Mij * source[j];
          }
          b[gsg](i) += temp;
        }

        // Solve system
        GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes));
      } // for gsg

      // Update phi
      auto& output_phi = GetDestinationPhi();
      for (int m = 0; m < num_moments_; ++m)
      {
        const double wn_d2m = d2m_op[m][direction_num];
        for (int i = 0; i < cell_num_nodes; ++i)
        {
          const size_t ir = cell_transport_view.MapDOF(i, m, gs_gi);
          for (int gsg = 0; gsg < gs_ss_size; ++gsg)
            output_phi[ir + gsg] += wn_d2m * b[gsg](i);
        }
      }

      // Save angular flux during sweep
      if (save_angular_flux_)
      {
        auto& output_psi = GetDestinationPsi();
        double* cell_psi_data =
          &output_psi[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];

        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const size_t imap =
            i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_ + gs_ss_begin;
          for (int gsg = 0; gsg < gs_ss_size; ++gsg)
            cell_psi_data[imap + gsg] = b[gsg](i);
        }
      }

      // For outoing, non-boundary faces, copy angular flux to fluds and
      // accumulate outflow
      int out_face_counter = -1;
      for (int f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::OUTGOING)
          continue;

        out_face_counter++;
        const auto& face = cell.faces[f];
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary_face = not face.has_neighbor;
        const bool is_reflecting_boundary_face =
          (is_boundary_face and angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting());
        const auto& IntF_shapeI = unit_cell_matrices_[cell_local_id].intS_shapeI[f];

        if (not is_boundary_face and not is_local_face)
          ++deploc_face_counter;

        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (int fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);

          if (is_boundary_face)
          {
            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              cell_transport_view.AddOutflow(
                f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg](i) * IntF_shapeI(i));
          }

          double* psi = nullptr;
          if (is_local_face)
            psi = fluds.OutgoingPsi(spls_index, out_face_counter, fi, as_ss_idx);
          else if (not is_boundary_face)
            psi = fluds.NLOutgoingPsi(deploc_face_counter, fi, as_ss_idx);
          else if (is_reflecting_boundary_face)
            psi = angle_set.PsiReflected(
              face.neighbor_id, direction_num, cell_local_id, f, fi, gs_ss_begin);
          else
            continue;

          if (not is_boundary_face or is_reflecting_boundary_face)
          {
            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              psi[gsg] = b[gsg](i);
          }
        } // for fi
      }   // for face
    }     // for angleset/subset
  }       // for cell
}

} // namespace opensn
