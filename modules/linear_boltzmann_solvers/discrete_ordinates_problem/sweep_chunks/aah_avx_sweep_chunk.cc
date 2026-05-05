// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_kernels.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"
#include <algorithm>
#include <array>
#include <vector>

namespace opensn
{

template <unsigned int NumNodes, bool time_dependent>
void
AAH_Sweep_FixedN(AAHSweepData& data, AngleSet& angle_set)
{
  static_assert(NumNodes >= 1 and NumNodes <= 8);

  CALI_CXX_MARK_SCOPE("AAH_Sweep_FixedN");

  const auto& groupset = data.groupset;
  const auto gs_size = groupset.GetNumGroups();
  const auto gs_gi = groupset.first_group;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();

  std::vector<double> b(static_cast<std::size_t>(gs_size) * NumNodes, 0.0);
  std::vector<double> sigma_block;
  sigma_block.reserve(data.group_block_size);

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const size_t num_spls = spls.size();
  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    const uint64_t cell_local_id = spls[spls_index];
    auto& cell = data.grid->local_cells[cell_local_id];
    auto& cell_transport_view = data.cell_transport_views[cell_local_id];
    const auto& cell_mapping = data.discretization.GetCellMapping(cell);
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();
    constexpr auto expected_nodes = static_cast<size_t>(NumNodes);
    OpenSnInvalidArgumentIf(cell_num_nodes != expected_nodes,
                            "AAH_Sweep_FixedN invoked for an incompatible cell topology.");

    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    const size_t cell_num_faces = cell.faces.size();
    std::vector<double> face_mu_values(cell_num_faces, 0.0);

    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;

    const auto& sigma_t = data.xs.at(cell.block_id)->GetSigmaTotal();

    const auto& unit_mats = data.unit_cell_matrices[cell_local_id];
    const auto& G = unit_mats.intV_shapeI_gradshapeJ;
    const auto& M = unit_mats.intV_shapeI_shapeJ;
    const auto& M_surf = unit_mats.intS_shapeI_shapeJ;

    const size_t matrix_size = static_cast<size_t>(NumNodes) * static_cast<size_t>(NumNodes);
    auto idx = [](int i, int j) -> size_t
    { return static_cast<size_t>(i) * static_cast<size_t>(NumNodes) + static_cast<size_t>(j); };

    std::array<double, matrix_size> mass_matrix{};
    PRAGMA_UNROLL
    for (int i = 0; i < NumNodes; ++i)
    {
      PRAGMA_UNROLL
      for (int j = 0; j < NumNodes; ++j)
        mass_matrix[idx(i, j)] = M(i, j);
    }

    std::array<double, matrix_size> Amat{};
    std::vector<std::array<size_t, NumNodes>> moment_dof_map(data.num_moments);
    for (unsigned int m = 0; m < data.num_moments; ++m)
    {
      PRAGMA_UNROLL
      for (int i = 0; i < NumNodes; ++i)
        moment_dof_map[m][i] = cell_transport_view.MapDOF(i, m, gs_gi);
    }

    std::vector<double> tau_gsg;
    if constexpr (time_dependent)
    {
      const auto& inv_velg = data.xs.at(cell.block_id)->GetInverseVelocity();
      const double theta = data.problem.GetTheta();
      const double inv_theta = 1.0 / theta;
      const double dt = data.problem.GetTimeStep();
      const double inv_dt = 1.0 / dt;
      tau_gsg.assign(gs_size, 0.0);
      for (size_t gsg = 0; gsg < gs_size; ++gsg)
        tau_gsg[gsg] = inv_velg[gs_gi + gsg] * inv_theta * inv_dt;
    }

    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      const auto direction_num = as_angle_indices[as_ss_idx];
      const auto omega = groupset.quadrature->omegas[direction_num];
      const auto wt = groupset.quadrature->weights[direction_num];

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      std::fill(b.begin(), b.end(), 0.0);

      PRAGMA_UNROLL
      for (int i = 0; i < NumNodes; ++i)
      {
        PRAGMA_UNROLL
        for (int j = 0; j < NumNodes; ++j)
          Amat[idx(i, j)] = omega.Dot(G(i, j));
      }

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

        const auto& Ms_f = M_surf[f];
        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        const double mu_f = -face_mu_values[f];

        for (size_t fj = 0; fj < num_face_nodes; ++fj)
        {
          const int j = cell_mapping.MapFaceNode(f, fj);
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
                                        data.surface_source_active);

          for (size_t fi = 0; fi < num_face_nodes; ++fi)
          {
            const int i = cell_mapping.MapFaceNode(f, fi);
            const double mu_Nij = mu_f * Ms_f(i, j);
            Amat[idx(i, j)] += mu_Nij;

            if (not psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg * NumNodes + i] += psi[gsg] * mu_Nij;
          }
        }
      }

      const auto dir_moment_offset =
        static_cast<std::size_t>(direction_num) * static_cast<std::size_t>(data.num_moments);
      const double* __restrict m2d_row = m2d_op.data() + dir_moment_offset;
      const double* __restrict d2m_row = d2m_op.data() + dir_moment_offset;

      const double* psi_old =
        (time_dependent and data.psi_old)
          ? &(*data.psi_old)[data.discretization.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)]
          : nullptr;

      for (unsigned int g0 = 0; g0 < gs_size; g0 += data.group_block_size)
      {
        const auto g1 = std::min(g0 + data.group_block_size, gs_size);
        const auto block_len = g1 - g0;
        sigma_block.resize(block_len);

        for (unsigned int gsg = g0; gsg < g1; ++gsg)
        {
          const size_t rel = gsg - g0;
          double sigma_tg = sigma_t[gs_gi + gsg];
          if constexpr (time_dependent)
            sigma_tg += tau_gsg[gsg];
          sigma_block[rel] = sigma_tg;

          double* __restrict bg = &b[static_cast<std::size_t>(gsg) * NumNodes];
          for (unsigned int m = 0; m < data.num_moments; ++m)
          {
            const double w = m2d_row[m];
            std::array<double, NumNodes> nodal_source{};
            for (int i = 0; i < NumNodes; ++i)
              nodal_source[i] = w * data.source_moments[moment_dof_map[m][i] + gsg];

            for (int i = 0; i < NumNodes; ++i)
            {
              double value = 0.0;
              const double* row = &mass_matrix[idx(i, 0)];
              PRAGMA_UNROLL
              for (int j = 0; j < NumNodes; ++j)
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

              for (int i = 0; i < NumNodes; ++i)
              {
                double value = 0.0;
                const double* row = &mass_matrix[idx(i, 0)];
                PRAGMA_UNROLL
                for (int j = 0; j < NumNodes; ++j)
                {
                  const size_t imap = static_cast<size_t>(j) * data.groupset_angle_group_stride +
                                      direction_num * data.groupset_group_stride;
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
          for (int i = 0; i < NumNodes; ++i)
          {
            PRAGMA_UNROLL
            for (int j = 0; j < NumNodes; ++j)
              A[idx(i, j)] = Amat[idx(i, j)] + sigma_tg * mass_matrix[idx(i, j)];
          }

          double* __restrict bg = &b[gsg * NumNodes];

          for (int pivot = 0; pivot < NumNodes; ++pivot)
          {
            const double inv = 1.0 / A[idx(pivot, pivot)];
            for (int row = pivot + 1; row < NumNodes; ++row)
            {
              const double factor = A[idx(row, pivot)] * inv;
              bg[row] -= factor * bg[pivot];
              PRAGMA_UNROLL
              for (int col = pivot + 1; col < NumNodes; ++col)
                A[idx(row, col)] -= factor * A[idx(pivot, col)];
            }
          }

          for (int pivot = NumNodes - 1; pivot >= 0; --pivot)
          {
            PRAGMA_UNROLL
            for (int col = pivot + 1; col < NumNodes; ++col)
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
            for (int i = 0; i < NumNodes; ++i)
            {
              const size_t dof = cell_transport_view.MapDOF(i, m, gs_gi);
              data.destination_phi[dof + gsg] += w * bg[i];
            }
          }
        }
      }

      if (data.save_angular_flux)
      {
        double* cell_psi_data =
          &data
             .destination_psi[data.discretization.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)];
        PRAGMA_UNROLL
        for (int i = 0; i < NumNodes; ++i)
        {
          const size_t imap =
            i * data.groupset_angle_group_stride + direction_num * data.groupset_group_stride;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
          {
            const double psi_sol = b[gsg * NumNodes + i];
            if constexpr (time_dependent)
            {
              const double theta = data.problem.GetTheta();
              const double inv_theta = 1.0 / theta;
              const double psi_old_val = psi_old ? psi_old[imap + gsg] : 0.0;
              cell_psi_data[imap + gsg] = inv_theta * (psi_sol + (theta - 1.0) * psi_old_val);
            }
            else
              cell_psi_data[imap + gsg] = psi_sol;
          }
        }
      }

      int out_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::OUTGOING)
          continue;

        const auto& face = cell.faces[f];
        const auto& IntF_shapeI = unit_mats.intS_shapeI[f];
        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary = not face.has_neighbor;
        const bool is_reflecting =
          is_boundary && angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting();
        const double mu_wt_f = wt * face_mu_values[f];

        ++out_face_counter;
        if (not is_boundary and not is_local_face)
          ++deploc_face_counter;

        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);

          if (is_boundary)
          {
            const double flux_i = mu_wt_f * IntF_shapeI(i);

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              cell_transport_view.AddOutflow(f, gs_gi + gsg, flux_i * b[gsg * NumNodes + i]);
          }

          double* psi = nullptr;
          if (is_local_face)
            psi = fluds.OutgoingPsi(spls_index, out_face_counter, fi, as_ss_idx);
          else if (not is_boundary)
            psi = fluds.NLOutgoingPsi(deploc_face_counter, fi, as_ss_idx);
          else if (is_reflecting)
            psi = angle_set.PsiReflected(face.neighbor_id, direction_num, cell_local_id, f, fi);
          else
            continue;

          if (not is_boundary or is_reflecting)
          {
            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              psi[gsg] = b[gsg * NumNodes + i];
          }
        }
      }
    }
  }
}

template <unsigned int NumNodes>
void
AAHSweepChunk::Sweep_FixedN(AngleSet& angle_set)
{
  AAHSweepData data{grid_,
                    discretization_,
                    unit_cell_matrices_,
                    cell_transport_views_,
                    source_moments_,
                    groupset_,
                    xs_,
                    num_moments_,
                    max_num_cell_dofs_,
                    min_num_cell_dofs_,
                    SaveAngularFluxEnabled(),
                    groupset_angle_group_stride_,
                    groupset_group_stride_,
                    destination_phi_,
                    destination_psi_,
                    surface_source_active_,
                    include_rhs_time_term_,
                    problem_,
                    nullptr,
                    group_block_size_};

  AAH_Sweep_FixedN<NumNodes, false>(data, angle_set);
}

template void AAH_Sweep_FixedN<2, false>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<3, false>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<4, false>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<5, false>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<6, false>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<7, false>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<8, false>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<2, true>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<3, true>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<4, true>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<5, true>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<6, true>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<7, true>(AAHSweepData&, AngleSet&);
template void AAH_Sweep_FixedN<8, true>(AAHSweepData&, AngleSet&);

template void AAHSweepChunk::Sweep_FixedN<2>(AngleSet&);
template void AAHSweepChunk::Sweep_FixedN<3>(AngleSet&);
template void AAHSweepChunk::Sweep_FixedN<4>(AngleSet&);
template void AAHSweepChunk::Sweep_FixedN<5>(AngleSet&);
template void AAHSweepChunk::Sweep_FixedN<6>(AngleSet&);
template void AAHSweepChunk::Sweep_FixedN<7>(AngleSet&);
template void AAHSweepChunk::Sweep_FixedN<8>(AngleSet&);

} // namespace opensn
