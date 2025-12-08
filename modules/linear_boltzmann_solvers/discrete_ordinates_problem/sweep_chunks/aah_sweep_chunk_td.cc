// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/logging/log_exceptions.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk_td.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"
#include <stdexcept>

namespace opensn
{

AAHSweepChunkTD::AAHSweepChunkTD(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
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
    problem_(problem),
    max_level_size_(problem.GetMaxLevelSize()),
    psi_old_(problem.GetPsiOldLocal()[groupset.id]),
    use_gpus_(problem.UseGPUs())
{
  if (use_gpus_)
    throw std::runtime_error("Time-dependent calculations do not yet support GPUs.\n");
}

void
AAHSweepChunkTD::Sweep(AngleSet& angle_set)
{
  CPUSweep(angle_set);
}

void
AAHSweepChunkTD::CPUSweep(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("AAHSweepChunkTD::Sweep");

  auto gs_size = groupset_.groups.size();
  auto gs_gi = groupset_.groups.front().id;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();

  DenseMatrix<double> Amat(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Atemp(max_num_cell_dofs_, max_num_cell_dofs_);
  std::vector<Vector<double>> b(groupset_.groups.size(), Vector<double>(max_num_cell_dofs_, 0.));
  std::vector<double> source(max_num_cell_dofs_);

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const size_t num_spls = spls.size();
  const auto& groups = groupset_.groups;

  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    auto cell_local_id = spls[spls_index];
    auto& cell = grid_->local_cells[cell_local_id];
    auto& cell_transport_view = cell_transport_views_[cell_local_id];
    const auto& cell_mapping = discretization_.GetCellMapping(cell);
    auto cell_num_faces = cell.faces.size();
    auto cell_num_nodes = cell_mapping.GetNumNodes();

    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    std::vector<double> face_mu_values(cell_num_faces);

    const auto& rho = densities_[cell.local_id];
    const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();

    const auto& inv_velg = xs_.at(cell.block_id)->GetInverseVelocity();
    const double theta = problem_.GetTheta();
    const double inv_theta = 1.0 / theta;
    const double dt = problem_.GetTimeStep();
    const double inv_dt = 1.0 / dt;
    std::vector<double> tau_gsg(gs_size, 0.0);
    for (size_t gsg = 0; gsg < gs_size; ++gsg)
      tau_gsg[gsg] = inv_velg[gs_gi + gsg] * inv_theta * inv_dt;

    const auto& G = unit_cell_matrices_[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = unit_cell_matrices_[cell_local_id].intS_shapeI_shapeJ;

    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;
    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      auto direction_num = as_angle_indices[as_ss_idx];
      auto omega = groupset_.quadrature->omegas[direction_num];
      auto wt = groupset_.quadrature->weights[direction_num];

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

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
                                          gs_gi,
                                          IsSurfaceSourceActive());

            if (not psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
          }
        }
      }

      const double* psi_old =
        &psi_old_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
      for (size_t gsg = 0; gsg < gs_size; ++gsg)
      {
        double sigma_tg = rho * sigma_t[gs_gi + gsg] + tau_gsg[gsg];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          double temp_src = 0.0;
          for (std::size_t m = 0; m < num_moments_; ++m)
          {
            const size_t ir = cell_transport_view.MapDOF(i, m, gs_gi + gsg);
            temp_src += m2d_op[direction_num][m] * source_moments_[ir];
          }
          const size_t imap =
            i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
          source[i] = temp_src;
          if (include_rhs_time_term_)
            source[i] += tau_gsg[gsg] * psi_old[imap + gsg];
        }

        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          double temp = 0.0;
          for (size_t j = 0; j < cell_num_nodes; ++j)
          {
            auto Mij = M(i, j);
            Atemp(i, j) = Amat(i, j) + Mij * sigma_tg;
            temp += Mij * source[j];
          }
          b[gsg](i) += temp;
        }

        GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes));
      }

      for (std::size_t m = 0; m < num_moments_; ++m)
      {
        const double wn_d2m = d2m_op[direction_num][m];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const size_t ir = cell_transport_view.MapDOF(i, m, gs_gi);
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            destination_phi_[ir + gsg] += wn_d2m * b[gsg](i);
        }
      }

      if (save_angular_flux_)
      {
        double* psi_new =
          &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const size_t imap =
            i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            psi_new[imap + gsg] = inv_theta * (b[gsg](i) + (theta - 1.0) * psi_old[imap + gsg]);
        }
      }

      int out_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
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

        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);

          if (is_boundary_face)
          {
            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              cell_transport_view.AddOutflow(
                f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg](i) * IntF_shapeI(i));
          }

          double* psi = nullptr;
          if (is_local_face)
            psi = fluds.OutgoingPsi(spls_index, out_face_counter, fi, as_ss_idx);
          else if (not is_boundary_face)
            psi = fluds.NLOutgoingPsi(deploc_face_counter, fi, as_ss_idx);
          else if (is_reflecting_boundary_face)
            psi = angle_set.PsiReflected(face.neighbor_id, direction_num, cell_local_id, f, fi);
          else
            continue;

          if (not is_boundary_face or is_reflecting_boundary_face)
          {
            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              psi[gsg] = b[gsg](i);
          }
        }
      }
    }
  }
}

AAHSweepChunkTD::~AAHSweepChunkTD() = default;

} // namespace opensn
