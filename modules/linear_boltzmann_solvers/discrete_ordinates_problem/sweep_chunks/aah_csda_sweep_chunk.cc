// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_csda_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/csda_sweep_helper.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "framework/data_types/dense_matrix.h"
#include "framework/data_types/vector.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"

namespace opensn
{

AAHCSDASweepChunk::AAHCSDASweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
  : SweepChunk(problem.GetPhiNewLocal(),
               problem.GetPsiNewLocal()[groupset.id],
               problem.GetGrid(),
               problem.GetSpatialDiscretization(),
               problem.GetUnitCellMatrices(),
               problem.GetCellTransportViews(),
               problem.GetCellOutflowViews(),
               problem.GetQMomentsLocal(),
               groupset,
               problem.GetBlockID2XSMap(),
               problem.GetNumMoments(),
               problem.GetMaxCellDOFCount(),
               problem.GetMinCellDOFCount()),
    problem_(problem),
    destination_phi_e_(problem.GetPhiENewLocal())
{
}

void
AAHCSDASweepChunk::ZeroDestinationPhi()
{
  SweepChunk::ZeroDestinationPhi();

  const auto gsi = groupset_.first_group;
  const auto gss = groupset_.GetNumGroups();
  const auto total_num_groups = problem_.GetNumGroups();

  for (const auto& cell : grid_->local_cells)
  {
    const auto mapping = cell.local_id * static_cast<size_t>(total_num_groups) + gsi;
    for (unsigned int g = 0; g < gss; ++g)
      destination_phi_e_[mapping + g] = 0.0;
  }
}

void
AAHCSDASweepChunk::Sweep(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("AAHCSDASweepChunk::Sweep");

  const auto gs_size = groupset_.GetNumGroups();
  const auto gs_gi = groupset_.first_group;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();

  DenseMatrix<double> Amat(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Atemp(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Aext(max_num_cell_dofs_ + 1, max_num_cell_dofs_ + 1);
  std::vector<Vector<double>> b(gs_size, Vector<double>(max_num_cell_dofs_, 0.0));
  Vector<double> b_ext(max_num_cell_dofs_ + 1, 0.0);
  std::vector<double> source(max_num_cell_dofs_);

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const size_t num_spls = spls.size();

  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    const auto cell_local_id = spls[spls_index];
    auto& cell = grid_->local_cells[cell_local_id];
    const auto& cell_transport_view = cell_transport_views_[cell_local_id];
    auto& cell_outflow_view = cell_outflow_views_[cell_local_id];
    const auto& cell_mapping = discretization_.GetCellMapping(cell);
    const size_t cell_num_faces = cell.faces.size();
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();

    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    std::vector<double> face_mu_values(cell_num_faces);

    const auto& xs = xs_.at(cell.block_id);
    const auto& sigma_t = xs->GetSigmaTotal();
    const auto csda_data = MakeCSDAMaterialData(*xs);

    const auto& unit_mats = unit_cell_matrices_[cell_local_id];
    const auto& G = unit_mats.intV_shapeI_gradshapeJ;
    const auto& M = unit_mats.intV_shapeI_shapeJ;
    const auto& v = unit_mats.intV_shapeI;
    const auto& M_surf = unit_mats.intS_shapeI_shapeJ;
    const auto& int_f_shape_i = unit_mats.intS_shapeI;

    double cell_volume = 0.0;
    for (size_t i = 0; i < cell_num_nodes; ++i)
      cell_volume += v(i);

    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;
    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();

    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      const auto direction_num = as_angle_indices[as_ss_idx];
      const auto omega = groupset_.quadrature->GetOmega(direction_num);
      const auto wt = groupset_.quadrature->GetWeight(direction_num);

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      std::vector<double> psiE_gsg(gs_size, 0.0);

      for (size_t gsg = 0; gsg < gs_size; ++gsg)
        for (size_t i = 0; i < cell_num_nodes; ++i)
          b[gsg](i) = 0.0;

      for (size_t i = 0; i < cell_num_nodes; ++i)
        for (size_t j = 0; j < cell_num_nodes; ++j)
          Amat(i, j) = omega.Dot(G(i, j));

      for (size_t f = 0; f < cell_num_faces; ++f)
        face_mu_values[f] = omega.Dot(cell.faces[f].normal);

      int in_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
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

        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          for (size_t fj = 0; fj < num_face_nodes; ++fj)
          {
            const int j = cell_mapping.MapFaceNode(f, fj);
            const double mu_Nij = -face_mu_values[f] * M_surf[f](i, j);
            Amat(i, j) += mu_Nij;

            const double* psi = nullptr;
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
                                          0,
                                          surface_source_active_);

            if (not psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
          }
        }
      }

      const double* m2d_row = m2d_op.data() + direction_num * static_cast<size_t>(num_moments_);
      const double* d2m_row = d2m_op.data() + direction_num * static_cast<size_t>(num_moments_);

      for (size_t gsg = 0; gsg < gs_size; ++gsg)
      {
        const double sigma_tg = sigma_t[gs_gi + gsg];

        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          double temp_src = 0.0;
          for (unsigned int m = 0; m < num_moments_; ++m)
          {
            const size_t ir = cell_transport_view.MapDOF(i, m, gs_gi + gsg);
            temp_src += m2d_row[m] * source_moments_[ir];
          }
          source[i] = temp_src;
        }

        if (IsCSDAActive(csda_data, gs_gi + gsg))
        {
          const CSDALocalSolveContext csda_ctx{gsg,
                                               gs_gi,
                                               sigma_tg,
                                               cell_volume,
                                               &Amat,
                                               &M,
                                               &v,
                                               &cell,
                                               &cell_mapping,
                                               &face_orientations,
                                               &face_mu_values,
                                               &int_f_shape_i,
                                               &cell_transport_view,
                                               &angle_set,
                                               ni_preloc_face_counter,
                                               &csda_data,
                                               &Aext,
                                               &b_ext,
                                               &psiE_gsg};
          SolveCSDALocalSystem(
            csda_ctx,
            [&source](size_t j) { return source[j]; },
            [&b](size_t gg, size_t j) { return b[gg](j); },
            [&b](size_t gg, size_t i, double value) { b[gg](i) = value; },
            [&fluds, &angle_set, &cell, cell_local_id, direction_num, spls_index, as_ss_idx](
              bool is_local,
              bool is_reflecting,
              size_t face_num,
              int in_face_counter_slope,
              int preloc_face_counter_slope,
              size_t gg)
            {
              const double* psiE = nullptr;
              if (is_local)
                psiE = fluds.UpwindPsiE(spls_index, in_face_counter_slope, 0, 0, as_ss_idx);
              else if (is_reflecting)
                psiE = angle_set.PsiBoundaryE(
                  cell.faces[face_num].neighbor_id, direction_num, cell_local_id, face_num, 0, 0);
              else
                psiE = fluds.NLUpwindPsiE(preloc_face_counter_slope, 0, 0, as_ss_idx);
              return psiE ? psiE[gg] : 0.0;
            });
        }
        else
        {
          for (size_t i = 0; i < cell_num_nodes; ++i)
          {
            double temp = 0.0;
            for (size_t j = 0; j < cell_num_nodes; ++j)
            {
              const double Mij = M(i, j);
              Atemp(i, j) = Amat(i, j) + Mij * sigma_tg;
              temp += Mij * source[j];
            }
            b[gsg](i) += temp;
          }

          GaussElimination(Atemp, b[gsg], cell_num_nodes);
        }
      }

      for (unsigned int m = 0; m < num_moments_; ++m)
      {
        const auto wn_d2m = d2m_row[m];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const size_t ir = cell_transport_view.MapDOF(i, m, gs_gi);
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            destination_phi_[ir + gsg] += wn_d2m * b[gsg](i);
        }
      }

      const double wn_d2m = d2m_row[0];
      const size_t cell_g_map =
        cell_local_id * static_cast<size_t>(problem_.GetNumGroups()) + gs_gi;
      for (size_t gsg = 0; gsg < gs_size; ++gsg)
        destination_phi_e_[cell_g_map + gsg] += wn_d2m * psiE_gsg[gsg];

      if (SaveAngularFluxEnabled())
      {
        double* psi_new =
          &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];

        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const size_t imap =
            i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            psi_new[imap + gsg] = b[gsg](i);
        }
      }

      int out_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::OUTGOING)
          continue;

        ++out_face_counter;
        const auto& face = cell.faces[f];
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary_face = not face.has_neighbor;
        const bool is_reflecting_boundary_face =
          (is_boundary_face and angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting());
        const auto& int_f_shape_i_face = int_f_shape_i[f];

        if (not is_boundary_face and not is_local_face)
          ++deploc_face_counter;

        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);

          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_outflow_view.Add(
              f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg](i) * int_f_shape_i_face(i));

          double* psi = nullptr;
          if (is_local_face)
            psi = fluds.OutgoingPsi(spls_index, out_face_counter, fi, as_ss_idx);
          else if (not is_boundary_face)
            psi = fluds.NLOutgoingPsi(deploc_face_counter, fi, as_ss_idx);
          else if (is_reflecting_boundary_face)
            psi = angle_set.PsiReflected(face.neighbor_id, direction_num, cell_local_id, f, fi);
          else
            continue;

          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            psi[gsg] = b[gsg](i);

          StoreOutgoingCSDASlope(
            (not is_boundary_face or is_reflecting_boundary_face),
            psiE_gsg,
            [&]() -> double*
            {
              if (is_local_face)
                return fluds.OutgoingPsiE(spls_index, out_face_counter, fi, as_ss_idx);
              if (is_reflecting_boundary_face)
                return angle_set.PsiReflectedE(
                  face.neighbor_id, direction_num, cell_local_id, f, fi);
              if (not is_boundary_face)
                return fluds.NLOutgoingPsiE(deploc_face_counter, fi, as_ss_idx);
              return nullptr;
            });
        }
      }
    }
  }
}

} // namespace opensn
