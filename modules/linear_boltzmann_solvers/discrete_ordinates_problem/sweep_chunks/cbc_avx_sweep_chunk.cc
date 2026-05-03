// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_kernels.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"
#include <algorithm>
#include <array>
#include <vector>

namespace opensn
{

template <unsigned int NumNodes, bool time_dependent>
void
CBC_Sweep_FixedN(CBCSweepData& data, AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBC_Sweep_FixedN");

  static_assert(NumNodes >= 2 and NumNodes <= 8);

  const auto& groupset = data.groupset;
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();

  OpenSnInvalidArgumentIf(data.cell_num_nodes != static_cast<size_t>(NumNodes),
                          "CBC_Sweep_FixedN invoked for an incompatible cell topology.");

  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[data.cell_local_id];
  const auto& sigma_t = data.xs.at(data.cell.block_id)->GetSigmaTotal();

  constexpr size_t matrix_size = static_cast<size_t>(NumNodes) * static_cast<size_t>(NumNodes);
  constexpr auto idx = [](size_t i, size_t j) -> size_t { return i * NumNodes + j; };

  std::array<double, matrix_size> mass_matrix{};
  PRAGMA_UNROLL
  for (size_t i = 0; i < NumNodes; ++i)
  {
    PRAGMA_UNROLL
    for (size_t j = 0; j < NumNodes; ++j)
      mass_matrix[idx(i, j)] = data.M(i, j);
  }

  auto& moment_dof_map = data.fixed_moment_dof_map;
  moment_dof_map.resize(static_cast<size_t>(data.num_moments) * NumNodes);
  for (unsigned int m = 0; m < data.num_moments; ++m)
  {
    PRAGMA_UNROLL
    for (size_t i = 0; i < NumNodes; ++i)
      moment_dof_map[static_cast<size_t>(m) * NumNodes + i] =
        data.cell_transport_view.MapDOF(i, m, data.gs_gi);
  }

  std::array<double, matrix_size> Amat{};
  auto& b = data.fixed_rhs_buffer;
  b.resize(static_cast<std::size_t>(data.gs_size) * NumNodes);
  auto& sigma_block = data.fixed_sigma_block;
  sigma_block.reserve(data.group_block_size);
  auto& face_mu_values = data.face_mu_values;
  face_mu_values.resize(data.cell_num_faces);

  std::vector<double> tau_gsg;
  if constexpr (time_dependent)
  {
    const auto& inv_velg = data.xs.at(data.cell.block_id)->GetInverseVelocity();
    const double theta = data.problem.GetTheta();
    const double inv_theta = 1.0 / theta;
    const double dt = data.problem.GetTimeStep();
    const double inv_dt = 1.0 / dt;

    tau_gsg.assign(data.gs_size, 0.0);
    for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
      tau_gsg[gsg] = inv_velg[data.gs_gi + gsg] * inv_theta * inv_dt;
  }

  const double* psi_old =
    (time_dependent and data.psi_old)
      ? &(*data.psi_old)[data.discretization.MapDOFLocal(data.cell, 0, groupset.psi_uk_man_, 0, 0)]
      : nullptr;

  const auto& as_angle_indices = angle_set.GetAngleIndices();
  auto& async_comm = data.async_comm;
  PrepareOutgoingNonlocalFaceBuffers(data, face_orientations);

  for (size_t as_ss_idx = 0; as_ss_idx < data.num_angles_in_as; ++as_ss_idx)
  {
    const auto direction_num = as_angle_indices[as_ss_idx];
    const auto& omega = groupset.quadrature->omegas[direction_num];
    const auto wt = groupset.quadrature->weights[direction_num];

    std::fill(b.begin(), b.end(), 0.0);

    PRAGMA_UNROLL
    for (size_t i = 0; i < NumNodes; ++i)
    {
      PRAGMA_UNROLL
      for (size_t j = 0; j < NumNodes; ++j)
        Amat[idx(i, j)] = omega.Dot(data.G(i, j));
    }

    for (size_t f = 0; f < data.cell_num_faces; ++f)
      face_mu_values[f] = omega.Dot(data.cell.faces[f].normal);

    for (size_t f = 0; f < data.cell_num_faces; ++f)
    {
      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const auto& face = data.cell.faces[f];
      const bool is_local_face = data.cell_transport_view.IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      const auto* face_nodal_mapping = is_boundary_face
                                         ? nullptr
                                         : &data.fluds.GetCommonData().GetFaceNodalMapping(
                                             data.cell_local_id, static_cast<unsigned int>(f));
      const auto incoming_nonlocal_slot =
        (is_boundary_face or is_local_face)
          ? CBC_FLUDSCommonData::INVALID_FACE_SLOT
          : data.fluds.GetCommonData().GetIncomingNonlocalFaceSlotByLocalFace(
              data.cell_local_id, static_cast<unsigned int>(f));

      const auto& Ms_f = data.M_surf[f];
      const size_t num_face_nodes = data.cell_mapping.GetNumFaceNodes(f);
      const double mu_f = -face_mu_values[f];

      for (size_t fj = 0; fj < num_face_nodes; ++fj)
      {
        const int j = data.cell_mapping.MapFaceNode(f, fj);

        const double* psi = nullptr;
        if (is_local_face)
          psi = data.fluds.UpwindPsi(*data.cell_transport_view.FaceNeighbor(f),
                                     face_nodal_mapping->cell_node_mapping_[fj],
                                     as_ss_idx);
        else if (not is_boundary_face)
          psi = data.fluds.NLUpwindPsi(
            incoming_nonlocal_slot, face_nodal_mapping->face_node_mapping_[fj], as_ss_idx);
        else
          psi = angle_set.PsiBoundary(face.neighbor_id,
                                      direction_num,
                                      data.cell_local_id,
                                      f,
                                      fj,
                                      data.gs_gi,
                                      data.surface_source_active);

        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = data.cell_mapping.MapFaceNode(f, fi);
          const double mu_Nij = mu_f * Ms_f(i, j);
          Amat[idx(i, j)] += mu_Nij;

          if (not psi)
            continue;

          for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
            b[gsg * NumNodes + i] += psi[gsg] * mu_Nij;
        }
      }
    }

    const auto dir_moment_offset =
      static_cast<std::size_t>(direction_num) * static_cast<std::size_t>(data.num_moments);
    const double* __restrict m2d_row = m2d_op.data() + dir_moment_offset;
    const double* __restrict d2m_row = d2m_op.data() + dir_moment_offset;

    for (unsigned int g0 = 0; g0 < data.gs_size; g0 += data.group_block_size)
    {
      const auto g1 = std::min(g0 + data.group_block_size, static_cast<unsigned int>(data.gs_size));
      const auto block_len = g1 - g0;
      sigma_block.resize(block_len);

      for (unsigned int gsg = g0; gsg < g1; ++gsg)
      {
        const size_t rel = gsg - g0;
        double sigma_tg = sigma_t[data.gs_gi + gsg];
        if constexpr (time_dependent)
          sigma_tg += tau_gsg[gsg];
        sigma_block[rel] = sigma_tg;

        auto* __restrict bg = &b[static_cast<std::size_t>(gsg) * NumNodes];
        std::array<double, NumNodes> nodal_source{};
        for (unsigned int m = 0; m < data.num_moments; ++m)
        {
          const double w = m2d_row[m];
          for (size_t i = 0; i < NumNodes; ++i)
            nodal_source[i] +=
              w * data.source_moments[moment_dof_map[static_cast<size_t>(m) * NumNodes + i] + gsg];
        }

        if constexpr (time_dependent)
        {
          if (data.include_rhs_time_term and psi_old)
          {
            const double tau = tau_gsg[gsg];
            for (size_t j = 0; j < NumNodes; ++j)
            {
              const size_t imap =
                j * data.groupset_angle_group_stride + direction_num * data.groupset_group_stride;
              nodal_source[j] += tau * psi_old[imap + gsg];
            }
          }
        }

        PRAGMA_UNROLL
        for (size_t i = 0; i < NumNodes; ++i)
        {
          double value = 0.0;
          const auto* row = &mass_matrix[idx(i, 0)];
          PRAGMA_UNROLL
          for (size_t j = 0; j < NumNodes; ++j)
            value += row[j] * nodal_source[j];
          bg[i] += value;
        }
      }

      size_t k = 0;

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
        const size_t gsg = g0 + k;
        const double sigma_tg = sigma_block[k];

        std::array<double, matrix_size> A{};
        PRAGMA_UNROLL
        for (size_t i = 0; i < NumNodes; ++i)
        {
          PRAGMA_UNROLL
          for (size_t j = 0; j < NumNodes; ++j)
            A[idx(i, j)] = Amat[idx(i, j)] + sigma_tg * mass_matrix[idx(i, j)];
        }

        auto* __restrict bg = &b[gsg * NumNodes];

        for (size_t pivot = 0; pivot < NumNodes; ++pivot)
        {
          const double inv = 1.0 / A[idx(pivot, pivot)];
          for (size_t row = pivot + 1; row < NumNodes; ++row)
          {
            const double factor = A[idx(row, pivot)] * inv;
            bg[row] -= factor * bg[pivot];
            PRAGMA_UNROLL
            for (size_t col = pivot + 1; col < NumNodes; ++col)
              A[idx(row, col)] -= factor * A[idx(pivot, col)];
          }
        }

        for (size_t pivot = NumNodes; pivot-- > 0;)
        {
          PRAGMA_UNROLL
          for (size_t col = pivot + 1; col < NumNodes; ++col)
            bg[pivot] -= A[idx(pivot, col)] * bg[col];
          bg[pivot] /= A[idx(pivot, pivot)];
        }
      }

      for (size_t gsg = g0; gsg < g1; ++gsg)
      {
        const auto* __restrict bg = &b[gsg * NumNodes];
        for (unsigned int m = 0; m < data.num_moments; ++m)
        {
          const double w = d2m_row[m];
          PRAGMA_UNROLL
          for (size_t i = 0; i < NumNodes; ++i)
          {
            const size_t dof = moment_dof_map[static_cast<size_t>(m) * NumNodes + i];
            data.destination_phi[dof + gsg] += w * bg[i];
          }
        }
      }
    }

    if (data.save_angular_flux)
    {
      auto* psi_new = &data.destination_psi[data.discretization.MapDOFLocal(
        data.cell, 0, groupset.psi_uk_man_, 0, 0)];

      double theta = 1.0;
      double inv_theta = 1.0;
      if constexpr (time_dependent)
      {
        theta = data.problem.GetTheta();
        inv_theta = 1.0 / theta;
      }

      PRAGMA_UNROLL
      for (size_t i = 0; i < NumNodes; ++i)
      {
        const size_t imap =
          i * data.groupset_angle_group_stride + direction_num * data.groupset_group_stride;

        for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
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

    for (size_t f = 0; f < data.cell_num_faces; ++f)
    {
      if (face_orientations[f] != FaceOrientation::OUTGOING)
        continue;

      const auto& face = data.cell.faces[f];
      const bool is_local_face = data.cell_transport_view.IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      const bool is_reflecting_boundary_face =
        (is_boundary_face and angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting());
      const auto& IntF_shapeI = data.IntS_shapeI[f];

      const size_t num_face_nodes = data.cell_mapping.GetNumFaceNodes(f);
      std::vector<double>* psi_nonlocal_outgoing = nullptr;
      if (not is_boundary_face and not is_local_face)
        psi_nonlocal_outgoing = &data.outgoing_nonlocal_face_buffer_by_face[f]->data;

      const double mu_wt_f = wt * face_mu_values[f];

      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = data.cell_mapping.MapFaceNode(f, fi);

        if (is_boundary_face)
        {
          const double flux_i = mu_wt_f * IntF_shapeI(i);
          for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
            data.cell_transport_view.AddOutflow(
              f, data.gs_gi + gsg, flux_i * b[gsg * NumNodes + i]);
        }

        double* psi = nullptr;
        if (is_local_face)
          psi = data.fluds.OutgoingPsi(data.cell, i, as_ss_idx);
        else if (not is_boundary_face)
          psi = data.fluds.NLOutgoingPsi(psi_nonlocal_outgoing, fi, as_ss_idx);
        else if (is_reflecting_boundary_face)
          psi = angle_set.PsiReflected(face.neighbor_id, direction_num, data.cell_local_id, f, fi);

        if (psi != nullptr)
        {
          for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
            psi[gsg] = b[gsg * NumNodes + i];
        }
      }
    }
  }

  QueueOutgoingNonlocalFaceBuffers(data, async_comm);
}

template <unsigned int NumNodes>
void
CBCSweepChunk::Sweep_FixedN(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunk::Sweep_FixedN");

  auto data = MakeCBCSweepData(discretization_,
                               source_moments_,
                               groupset_,
                               xs_,
                               num_moments_,
                               max_num_cell_dofs_,
                               SaveAngularFluxEnabled(),
                               groupset_angle_group_stride_,
                               groupset_group_stride_,
                               destination_phi_,
                               destination_psi_,
                               include_rhs_time_term_,
                               problem_,
                               nullptr,
                               group_block_size_,
                               ctx_);

  CBC_Sweep_FixedN<NumNodes, false>(data, angle_set);
}

template void CBC_Sweep_FixedN<2, false>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<3, false>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<4, false>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<5, false>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<6, false>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<7, false>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<8, false>(CBCSweepData&, AngleSet&);

template void CBC_Sweep_FixedN<2, true>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<3, true>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<4, true>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<5, true>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<6, true>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<7, true>(CBCSweepData&, AngleSet&);
template void CBC_Sweep_FixedN<8, true>(CBCSweepData&, AngleSet&);

template void CBCSweepChunk::Sweep_FixedN<2>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<3>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<4>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<5>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<6>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<7>(AngleSet&);
template void CBCSweepChunk::Sweep_FixedN<8>(AngleSet&);

} // namespace opensn
