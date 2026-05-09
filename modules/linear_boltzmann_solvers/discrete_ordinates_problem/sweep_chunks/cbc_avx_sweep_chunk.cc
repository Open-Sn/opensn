// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk_td.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_kernels.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "caliper/cali.h"
#include <algorithm>
#include <array>
#include <vector>

namespace opensn
{

template <unsigned int NumNodes, bool time_dependent, class SweepChunkT>
void
CBC_Sweep_FixedN(SweepChunkT& sweep_chunk, AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBC_Sweep_FixedN");

  static_assert(NumNodes >= 2 and NumNodes <= 8);

  const auto& groupset = sweep_chunk.groupset_;
  const auto gs_size = groupset.GetNumGroups();
  const auto gs_gi = groupset.first_group;
  const auto num_angles_in_as = angle_set.GetNumAngles();
  const auto group_angle_stride = gs_size * num_angles_in_as;
  const auto& cell = *sweep_chunk.cell_;
  const auto cell_local_id = cell.local_id;
  const auto& cell_mapping = sweep_chunk.discretization_.GetCellMapping(cell);
  const auto& cell_transport_view = sweep_chunk.cell_transport_views_[cell_local_id];
  auto& cell_outflow_view = sweep_chunk.cell_outflow_views_[cell_local_id];
  const std::size_t cell_num_faces = cell.faces.size();
  const auto& unit_mats = sweep_chunk.unit_cell_matrices_[cell_local_id];
  auto& fluds = *sweep_chunk.fluds_;
  const auto& common_data = fluds.GetCommonData();
  const auto& G = unit_mats.intV_shapeI_gradshapeJ;
  const auto& M = unit_mats.intV_shapeI_shapeJ;
  const auto& M_surf = unit_mats.intS_shapeI_shapeJ;
  const auto& IntS_shapeI = unit_mats.intS_shapeI;
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();

  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[cell_local_id];
  const auto& sigma_t = sweep_chunk.xs_.at(cell.block_id)->GetSigmaTotal();

  constexpr std::size_t matrix_size =
    static_cast<std::size_t>(NumNodes) * static_cast<std::size_t>(NumNodes);
  constexpr auto idx = [](std::size_t i, std::size_t j) -> std::size_t { return i * NumNodes + j; };

  std::array<double, matrix_size> mass_matrix{};
  PRAGMA_UNROLL
  for (std::size_t i = 0; i < NumNodes; ++i)
  {
    PRAGMA_UNROLL
    for (std::size_t j = 0; j < NumNodes; ++j)
      mass_matrix[idx(i, j)] = M(i, j);
  }

  auto& moment_dof_map = sweep_chunk.workspace_.moment_dof_map;
  moment_dof_map.resize(static_cast<std::size_t>(sweep_chunk.num_moments_) * NumNodes);
  for (unsigned int m = 0; m < sweep_chunk.num_moments_; ++m)
  {
    PRAGMA_UNROLL
    for (std::size_t i = 0; i < NumNodes; ++i)
      moment_dof_map[static_cast<std::size_t>(m) * NumNodes + i] =
        cell_transport_view.MapDOF(i, m, gs_gi);
  }

  std::array<double, matrix_size> Amat{};
  auto& b = sweep_chunk.workspace_.rhs_values;
  b.resize(static_cast<std::size_t>(gs_size) * NumNodes);
  auto& sigma_block = sweep_chunk.workspace_.sigma_t_group_block;
  sigma_block.reserve(sweep_chunk.group_block_size_);
  auto& face_mu_values = sweep_chunk.workspace_.face_mu_values;
  face_mu_values.resize(cell_num_faces);

  std::vector<double> tau_gsg;
  if constexpr (time_dependent)
  {
    const auto& inv_velg = sweep_chunk.xs_.at(cell.block_id)->GetInverseVelocity();
    const double theta = sweep_chunk.problem_.GetTheta();
    const double inv_theta = 1.0 / theta;
    const double dt = sweep_chunk.problem_.GetTimeStep();
    const double inv_dt = 1.0 / dt;

    tau_gsg.assign(gs_size, 0.0);
    for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
      tau_gsg[gsg] = inv_velg[gs_gi + gsg] * inv_theta * inv_dt;
  }

  const double* psi_old = nullptr;
  if constexpr (time_dependent)
    psi_old =
      &sweep_chunk
         .psi_old_[sweep_chunk.discretization_.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)];

  const auto& as_angle_indices = angle_set.GetAngleIndices();
  PrepareNonlocalOutgoingPsi(sweep_chunk.workspace_,
                             fluds,
                             cell,
                             cell_local_id,
                             cell_mapping,
                             cell_transport_view,
                             group_angle_stride,
                             face_orientations);

  for (std::size_t as_ss_idx = 0; as_ss_idx < num_angles_in_as; ++as_ss_idx)
  {
    const auto direction_num = as_angle_indices[as_ss_idx];
    const auto& omega = groupset.quadrature->GetOmega(direction_num);
    const auto wt = groupset.quadrature->GetWeight(direction_num);

    std::fill(b.begin(), b.end(), 0.0);

    PRAGMA_UNROLL
    for (std::size_t i = 0; i < NumNodes; ++i)
    {
      PRAGMA_UNROLL
      for (std::size_t j = 0; j < NumNodes; ++j)
        Amat[idx(i, j)] = omega.Dot(G(i, j));
    }

    for (std::size_t f = 0; f < cell_num_faces; ++f)
      face_mu_values[f] = omega.Dot(cell.faces[f].normal);

    for (std::size_t f = 0; f < cell_num_faces; ++f)
    {
      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const auto& face = cell.faces[f];
      const bool is_local_face = cell_transport_view.IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      const auto* face_nodal_mapping =
        is_boundary_face
          ? nullptr
          : &common_data.GetFaceNodalMapping(cell_local_id, static_cast<unsigned int>(f));
      const bool is_delayed_local_face =
        is_local_face and
        common_data.IsDelayedLocalIncomingFace(cell_local_id, static_cast<unsigned int>(f));
      const bool is_delayed_nonlocal_face =
        (not is_boundary_face) and (not is_local_face) and
        common_data.IsDelayedNonlocalIncomingFace(cell_local_id, static_cast<unsigned int>(f));
      const auto delayed_nonlocal_face_info =
        is_delayed_nonlocal_face
          ? common_data.DelayedNonlocalFaceByLocalFace(cell_local_id, static_cast<unsigned int>(f))
          : CBC_FLUDSCommonData::DelayedNonlocalFaceInfo{};
      const auto incoming_nonlocal_slot =
        (is_boundary_face or is_local_face or is_delayed_nonlocal_face)
          ? CBC_FLUDSCommonData::INVALID_FACE_SLOT
          : common_data.IncomingFaceSlot(cell_local_id, static_cast<unsigned int>(f));

      const auto& Ms_f = M_surf[f];
      const std::size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      const double mu_f = -face_mu_values[f];

      for (std::size_t fj = 0; fj < num_face_nodes; ++fj)
      {
        const int j = cell_mapping.MapFaceNode(f, fj);

        const double* psi = nullptr;
        if (is_delayed_local_face)
          psi = fluds.DelayedUpwindPsi(cell_local_id,
                                       static_cast<unsigned int>(f),
                                       face_nodal_mapping->face_node_mapping_[fj],
                                       as_ss_idx);
        else if (is_delayed_nonlocal_face)
          psi = fluds.DelayedNLUpwindPsi(
            delayed_nonlocal_face_info, face_nodal_mapping->face_node_mapping_[fj], as_ss_idx);
        else if (is_local_face)
          psi = fluds.UpwindPsi(*cell_transport_view.FaceNeighbor(f),
                                face_nodal_mapping->cell_node_mapping_[fj],
                                as_ss_idx);
        else if (not is_boundary_face)
          psi = fluds.NLUpwindPsi(
            incoming_nonlocal_slot, face_nodal_mapping->face_node_mapping_[fj], as_ss_idx);
        else
          psi = angle_set.PsiBoundary(face.neighbor_id,
                                      direction_num,
                                      cell_local_id,
                                      f,
                                      fj,
                                      0,
                                      sweep_chunk.IsSurfaceSourceActive());

        for (std::size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          const double mu_Nij = mu_f * Ms_f(i, j);
          Amat[idx(i, j)] += mu_Nij;

          if (not psi)
            continue;

          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            b[gsg * NumNodes + i] += psi[gsg] * mu_Nij;
        }
      }
    }

    const auto dir_moment_offset =
      static_cast<std::size_t>(direction_num) * static_cast<std::size_t>(sweep_chunk.num_moments_);
    const double* __restrict m2d_row = m2d_op.data() + dir_moment_offset;
    const double* __restrict d2m_row = d2m_op.data() + dir_moment_offset;

    for (unsigned int g0 = 0; g0 < gs_size; g0 += sweep_chunk.group_block_size_)
    {
      const auto g1 = std::min(g0 + sweep_chunk.group_block_size_, gs_size);
      const auto block_len = g1 - g0;
      sigma_block.resize(block_len);

      for (unsigned int gsg = g0; gsg < g1; ++gsg)
      {
        const std::size_t rel = gsg - g0;
        double sigma_tg = sigma_t[gs_gi + gsg];
        if constexpr (time_dependent)
          sigma_tg += tau_gsg[gsg];
        sigma_block[rel] = sigma_tg;

        auto* __restrict bg = &b[static_cast<std::size_t>(gsg) * NumNodes];
        std::array<double, NumNodes> nodal_source{};
        for (unsigned int m = 0; m < sweep_chunk.num_moments_; ++m)
        {
          const double w = m2d_row[m];
          for (std::size_t i = 0; i < NumNodes; ++i)
            nodal_source[i] +=
              w *
              sweep_chunk
                .source_moments_[moment_dof_map[static_cast<std::size_t>(m) * NumNodes + i] + gsg];
        }

        if constexpr (time_dependent)
        {
          if (sweep_chunk.include_rhs_time_term_ and psi_old)
          {
            const double tau = tau_gsg[gsg];
            for (std::size_t j = 0; j < NumNodes; ++j)
            {
              const std::size_t imap = j * sweep_chunk.groupset_angle_group_stride_ +
                                       direction_num * sweep_chunk.groupset_group_stride_;
              nodal_source[j] += tau * psi_old[imap + gsg];
            }
          }
        }

        PRAGMA_UNROLL
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
          double value = 0.0;
          const auto* row = &mass_matrix[idx(i, 0)];
          PRAGMA_UNROLL
          for (std::size_t j = 0; j < NumNodes; ++j)
            value += row[j] * nodal_source[j];
          bg[i] += value;
        }
      }

      std::size_t k = 0;

#if __AVX512F__
      for (; k + simd_width <= block_len; k += simd_width)
        detail::SimdBatchSolve<detail::AVX512Ops, NumNodes>(
          Amat.data(), mass_matrix.data(), &sigma_block[k], &b[(g0 + k) * NumNodes]);
#elif __AVX2__
      for (; k + simd_width <= block_len; k += simd_width)
        detail::SimdBatchSolve<detail::AVX2Ops, NumNodes>(
          Amat.data(), mass_matrix.data(), &sigma_block[k], &b[(g0 + k) * NumNodes]);
#endif

      for (; k < block_len; ++k)
      {
        const std::size_t gsg = g0 + k;
        const double sigma_tg = sigma_block[k];

        std::array<double, matrix_size> A{};
        PRAGMA_UNROLL
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
          PRAGMA_UNROLL
          for (std::size_t j = 0; j < NumNodes; ++j)
            A[idx(i, j)] = Amat[idx(i, j)] + sigma_tg * mass_matrix[idx(i, j)];
        }

        auto* __restrict bg = &b[gsg * NumNodes];

        for (std::size_t pivot = 0; pivot < NumNodes; ++pivot)
        {
          const double inv = 1.0 / A[idx(pivot, pivot)];
          for (std::size_t row = pivot + 1; row < NumNodes; ++row)
          {
            const double factor = A[idx(row, pivot)] * inv;
            bg[row] -= factor * bg[pivot];
            PRAGMA_UNROLL
            for (std::size_t col = pivot + 1; col < NumNodes; ++col)
              A[idx(row, col)] -= factor * A[idx(pivot, col)];
          }
        }

        for (std::size_t pivot = NumNodes; pivot-- > 0;)
        {
          PRAGMA_UNROLL
          for (std::size_t col = pivot + 1; col < NumNodes; ++col)
            bg[pivot] -= A[idx(pivot, col)] * bg[col];
          bg[pivot] /= A[idx(pivot, pivot)];
        }
      }

      for (std::size_t gsg = g0; gsg < g1; ++gsg)
      {
        const auto* __restrict bg = &b[gsg * NumNodes];
        for (unsigned int m = 0; m < sweep_chunk.num_moments_; ++m)
        {
          const double w = d2m_row[m];
          PRAGMA_UNROLL
          for (std::size_t i = 0; i < NumNodes; ++i)
          {
            const std::size_t dof = moment_dof_map[static_cast<std::size_t>(m) * NumNodes + i];
            sweep_chunk.destination_phi_[dof + gsg] += w * bg[i];
          }
        }
      }
    }

    if (sweep_chunk.SaveAngularFluxEnabled())
    {
      auto* psi_new = &sweep_chunk.destination_psi_[sweep_chunk.discretization_.MapDOFLocal(
        cell, 0, groupset.psi_uk_man_, 0, 0)];

      double theta = 1.0;
      double inv_theta = 1.0;
      if constexpr (time_dependent)
      {
        theta = sweep_chunk.problem_.GetTheta();
        inv_theta = 1.0 / theta;
      }

      PRAGMA_UNROLL
      for (std::size_t i = 0; i < NumNodes; ++i)
      {
        const std::size_t imap = i * sweep_chunk.groupset_angle_group_stride_ +
                                 direction_num * sweep_chunk.groupset_group_stride_;

        for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
        {
          const double psi_sol = b[gsg * NumNodes + i];
          if constexpr (time_dependent)
          {
            const double psi_old_val = psi_old ? psi_old[imap + gsg] : 0.0;
            psi_new[imap + gsg] = inv_theta * (psi_sol + (theta - 1.0) * psi_old_val);
          }
          else
            psi_new[imap + gsg] = psi_sol;
        }
      }
    }

    for (std::size_t f = 0; f < cell_num_faces; ++f)
    {
      if (face_orientations[f] != FaceOrientation::OUTGOING)
        continue;

      const auto& face = cell.faces[f];
      const bool is_local_face = cell_transport_view.IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      const bool is_reflecting_boundary_face =
        (is_boundary_face and angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting());
      const auto* face_nodal_mapping =
        is_boundary_face
          ? nullptr
          : &common_data.GetFaceNodalMapping(cell_local_id, static_cast<unsigned int>(f));
      const bool is_delayed_local_outgoing =
        is_local_face and
        common_data.IsDelayedLocalOutgoingFace(cell_local_id, static_cast<unsigned int>(f));
      std::uint32_t delayed_local_cell_local_id = 0;
      unsigned int delayed_local_face_id = 0;
      if (is_delayed_local_outgoing)
      {
        delayed_local_cell_local_id = face.GetNeighborLocalID(fluds.GetSPDS().GetGrid().get());
        delayed_local_face_id = static_cast<unsigned int>(face_nodal_mapping->associated_face_);
      }
      const auto& IntF_shapeI = IntS_shapeI[f];

      const std::size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      std::vector<double>* psi_nonlocal_outgoing = nullptr;
      if (not is_boundary_face and not is_local_face)
        psi_nonlocal_outgoing = &sweep_chunk.workspace_.outgoing_psi_by_face[f]->data;

      const double mu_wt_f = wt * face_mu_values[f];

      for (std::size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping.MapFaceNode(f, fi);

        if (is_boundary_face)
        {
          const double flux_i = mu_wt_f * IntF_shapeI(i);
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_outflow_view.Add(f, gs_gi + gsg, flux_i * b[gsg * NumNodes + i]);
        }

        double* psi = nullptr;
        if (is_delayed_local_outgoing)
          psi = fluds.DelayedLocalOutgoingPsi(delayed_local_cell_local_id,
                                              delayed_local_face_id,
                                              static_cast<unsigned int>(fi),
                                              as_ss_idx);
        else if (is_local_face)
          psi = fluds.OutgoingPsi(cell, i, as_ss_idx);
        else if (not is_boundary_face)
          psi = fluds.NLOutgoingPsi(psi_nonlocal_outgoing, fi, as_ss_idx);
        else if (is_reflecting_boundary_face)
          psi = angle_set.PsiReflected(face.neighbor_id, direction_num, cell_local_id, f, fi);

        if (psi != nullptr)
        {
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            psi[gsg] = b[gsg * NumNodes + i];
        }
      }
    }
  }

  QueueNonlocalOutgoingPsi(sweep_chunk.workspace_, *sweep_chunk.async_comm_);
}

template <unsigned int NumNodes>
void
CBCSweepChunk::Sweep_FixedN(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunk::Sweep_FixedN");

  CBC_Sweep_FixedN<NumNodes, false>(*this, angle_set);
}

template void CBC_Sweep_FixedN<2, false>(CBCSweepChunk&, AngleSet&);
template void CBC_Sweep_FixedN<3, false>(CBCSweepChunk&, AngleSet&);
template void CBC_Sweep_FixedN<4, false>(CBCSweepChunk&, AngleSet&);
template void CBC_Sweep_FixedN<5, false>(CBCSweepChunk&, AngleSet&);
template void CBC_Sweep_FixedN<6, false>(CBCSweepChunk&, AngleSet&);
template void CBC_Sweep_FixedN<7, false>(CBCSweepChunk&, AngleSet&);
template void CBC_Sweep_FixedN<8, false>(CBCSweepChunk&, AngleSet&);

template void CBC_Sweep_FixedN<2, true>(CBCSweepChunkTD&, AngleSet&);
template void CBC_Sweep_FixedN<3, true>(CBCSweepChunkTD&, AngleSet&);
template void CBC_Sweep_FixedN<4, true>(CBCSweepChunkTD&, AngleSet&);
template void CBC_Sweep_FixedN<5, true>(CBCSweepChunkTD&, AngleSet&);
template void CBC_Sweep_FixedN<6, true>(CBCSweepChunkTD&, AngleSet&);
template void CBC_Sweep_FixedN<7, true>(CBCSweepChunkTD&, AngleSet&);
template void CBC_Sweep_FixedN<8, true>(CBCSweepChunkTD&, AngleSet&);

template void CBCSweepChunk::Sweep_FixedN<2>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<3>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<4>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<5>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<6>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<7>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<8>(AngleSet&);

} // namespace opensn
