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

  struct IncomingFaceData
  {
    const FaceNodalMapping* face_nodal_mapping = nullptr;
    double* psi_base = nullptr;
  };

  struct OutgoingFaceData
  {
    bool is_reflecting_boundary_face = false;
    double* psi_base = nullptr;
    const CBC_FLUDSCommonData::OutgoingNonlocalFaceInfo* outgoing_nonlocal_face_info = nullptr;
  };

  const auto& groupset = data.groupset;
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();

  OpenSnInvalidArgumentIf(data.cell_num_nodes != static_cast<size_t>(NumNodes),
                          "CBC_Sweep_FixedN invoked for an incompatible cell topology.");

  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[data.cell_local_id];
  const auto& cell_xs = data.cell_transport_view.GetXS();
  const auto& sigma_t = cell_xs.GetSigmaTotal();

  constexpr size_t matrix_size = static_cast<size_t>(NumNodes) * static_cast<size_t>(NumNodes);
  auto idx = [](size_t i, size_t j) -> size_t { return i * NumNodes + j; };

  std::array<double, matrix_size> mass_matrix{};
  PRAGMA_UNROLL
  for (size_t i = 0; i < NumNodes; ++i)
  {
    PRAGMA_UNROLL
    for (size_t j = 0; j < NumNodes; ++j)
      mass_matrix[idx(i, j)] = data.M(i, j);
  }

  std::vector<std::array<size_t, NumNodes>> moment_dof_map(data.num_moments);
  for (unsigned int m = 0; m < data.num_moments; ++m)
  {
    PRAGMA_UNROLL
    for (size_t i = 0; i < NumNodes; ++i)
      moment_dof_map[m][i] = data.cell_transport_view.MapDOF(i, m, data.gs_gi);
  }

  std::array<double, matrix_size> Amat{};
  std::vector<double> b(static_cast<std::size_t>(data.gs_size) * NumNodes, 0.0);
  std::vector<double> sigma_block;
  sigma_block.reserve(data.group_block_size);
  std::vector<double> face_mu_values(data.cell_num_faces);

  std::vector<double> tau_gsg;
  if constexpr (time_dependent)
  {
    const auto& inv_velg = cell_xs.GetInverseVelocity();
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
  const auto& cbc_common = dynamic_cast<const CBC_FLUDSCommonData&>(data.fluds.GetCommonData());
  auto* const async_comm = dynamic_cast<CBC_AsynchronousCommunicator*>(angle_set.GetCommunicator());
  std::vector<IncomingFaceData> incoming_face_data(data.cell_num_faces);
  std::vector<OutgoingFaceData> outgoing_face_data(data.cell_num_faces);
  for (size_t f = 0; f < data.cell_num_faces; ++f)
  {
    const auto& face = data.cell.faces[f];
    const bool is_local_face = data.cell_transport_view.IsFaceLocal(f);
    const bool is_boundary_face = not face.has_neighbor;
    const auto* face_nodal_mapping =
      &data.fluds.GetCommonData().GetFaceNodalMapping(data.cell_local_id, f);

    if (face_orientations[f] == FaceOrientation::INCOMING)
    {
      auto& face_data = incoming_face_data[f];
      face_data.face_nodal_mapping = face_nodal_mapping;
      if (is_local_face)
        face_data.psi_base =
          data.fluds.GetLocalFacePsiPointer(data.cell_local_id, static_cast<unsigned int>(f));
      else if (not is_boundary_face)
        face_data.psi_base = data.fluds.GetIncomingNonlocalFacePsiPointer(
          data.cell_local_id, static_cast<unsigned int>(f));
    }

    if (face_orientations[f] == FaceOrientation::OUTGOING)
    {
      auto& face_data = outgoing_face_data[f];
      face_data.is_reflecting_boundary_face =
        is_boundary_face and angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting();
      if (is_local_face)
        face_data.psi_base =
          data.fluds.GetLocalFacePsiPointer(data.cell_local_id, static_cast<unsigned int>(f));
      if (not is_local_face and not is_boundary_face)
        face_data.outgoing_nonlocal_face_info =
          &cbc_common.GetOutgoingNonlocalFaceInfo(data.cell_local_id, static_cast<unsigned int>(f));
    }
  }

  double* psi_new_base = nullptr;
  double theta = 1.0;
  double inv_theta = 1.0;
  if (data.save_angular_flux)
  {
    psi_new_base = &data.destination_psi[data.discretization.MapDOFLocal(
      data.cell, 0, groupset.psi_uk_man_, 0, 0)];
    if constexpr (time_dependent)
    {
      theta = data.problem.GetTheta();
      inv_theta = 1.0 / theta;
    }
  }

  for (size_t as_ss_idx = 0; as_ss_idx < data.num_angles_in_as; ++as_ss_idx)
  {
    const auto direction_num = as_angle_indices[as_ss_idx];
    const auto omega = groupset.quadrature->omegas[direction_num];
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
      const auto& face_data = incoming_face_data[f];
      const auto* face_nodal_mapping = face_data.face_nodal_mapping;

      const auto& Ms_f = data.M_surf[f];
      const size_t num_face_nodes = data.cell_mapping.GetNumFaceNodes(f);
      const double mu_f = -face_mu_values[f];

      for (size_t fj = 0; fj < num_face_nodes; ++fj)
      {
        const int j = data.cell_mapping.MapFaceNode(f, fj);

        const double* psi = nullptr;
        if (face_data.psi_base != nullptr)
          psi = face_data.psi_base +
                static_cast<size_t>(face_nodal_mapping->face_node_mapping_[fj]) *
                  data.group_angle_stride +
                as_ss_idx * data.group_stride;
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

        double* __restrict bg = &b[static_cast<std::size_t>(gsg) * NumNodes];
        for (unsigned int m = 0; m < data.num_moments; ++m)
        {
          const double w = m2d_row[m];
          std::array<double, NumNodes> nodal_source{};
          for (size_t i = 0; i < NumNodes; ++i)
            nodal_source[i] = w * data.source_moments[moment_dof_map[m][i] + gsg];

          for (size_t i = 0; i < NumNodes; ++i)
          {
            double value = 0.0;
            const double* row = &mass_matrix[idx(i, 0)];
            PRAGMA_UNROLL
            for (size_t j = 0; j < NumNodes; ++j)
              value += row[j] * nodal_source[j];
            bg[i] += value;
          }
        }
      }

      if constexpr (time_dependent)
      {
        if (data.include_rhs_time_term and psi_old)
        {
          for (size_t gsg = g0; gsg < g1; ++gsg)
          {
            const double tau = tau_gsg[gsg];
            double* __restrict bg = &b[gsg * NumNodes];

            for (size_t i = 0; i < NumNodes; ++i)
            {
              double value = 0.0;
              const double* row = &mass_matrix[idx(i, 0)];
              PRAGMA_UNROLL
              for (size_t j = 0; j < NumNodes; ++j)
              {
                const size_t imap =
                  j * data.groupset_angle_group_stride + direction_num * data.groupset_group_stride;
                const double psi_old_val = psi_old[imap + gsg];
                value += row[j] * psi_old_val;
              }
              bg[i] += tau * value;
            }
          }
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

        double* __restrict bg = &b[gsg * NumNodes];

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
        const double* __restrict bg = &b[gsg * NumNodes];
        for (unsigned int m = 0; m < data.num_moments; ++m)
        {
          const double w = d2m_row[m];
          PRAGMA_UNROLL
          for (size_t i = 0; i < NumNodes; ++i)
            data.destination_phi[moment_dof_map[m][i] + gsg] += w * bg[i];
        }
      }
    }

    if (data.save_angular_flux)
    {
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
            psi_new_base[imap + gsg] = inv_theta * (psi_sol + (theta - 1.0) * psi_old_val);
          }
          else
            psi_new_base[imap + gsg] = psi_sol;
        }
      }
    }

    for (size_t f = 0; f < data.cell_num_faces; ++f)
    {
      if (face_orientations[f] != FaceOrientation::OUTGOING)
        continue;

      const auto& face = data.cell.faces[f];
      const auto& face_data = outgoing_face_data[f];
      const bool is_reflecting_boundary_face = face_data.is_reflecting_boundary_face;
      const auto& IntF_shapeI = data.IntS_shapeI[f];

      const size_t num_face_nodes = data.cell_mapping.GetNumFaceNodes(f);
      double* psi_nonlocal_outgoing = nullptr;

      if (face_data.outgoing_nonlocal_face_info != nullptr)
      {
        const auto& outgoing_nonlocal_face_info = *face_data.outgoing_nonlocal_face_info;
        const size_t data_size_for_msg =
          static_cast<size_t>(outgoing_nonlocal_face_info.num_face_nodes) * data.group_angle_stride;
        psi_nonlocal_outgoing =
          async_comm
            ->InitGetDownwindMessageData(outgoing_nonlocal_face_info.locality,
                                         outgoing_nonlocal_face_info.cell_global_id,
                                         outgoing_nonlocal_face_info.associated_face,
                                         angle_set.GetID(),
                                         data_size_for_msg)
            .data();
      }

      const double mu_wt_f = wt * face_mu_values[f];

      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = data.cell_mapping.MapFaceNode(f, fi);

        if (face_data.outgoing_nonlocal_face_info == nullptr and face_data.psi_base == nullptr)
        {
          const double flux_i = mu_wt_f * IntF_shapeI(i);
          for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
            data.cell_transport_view.AddOutflow(
              f, data.gs_gi + gsg, flux_i * b[gsg * NumNodes + i]);
        }

        double* psi = nullptr;
        if (face_data.psi_base != nullptr)
          psi = face_data.psi_base + fi * data.group_angle_stride + as_ss_idx * data.group_stride;
        else if (face_data.outgoing_nonlocal_face_info != nullptr)
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
}

template <unsigned int NumNodes>
void
CBCSweepChunk::Sweep_FixedN(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunk::Sweep_FixedN");

  auto data = MakeSweepData(nullptr);

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
