// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/dense_matrix.h"
#include "framework/data_types/vector.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_view.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include <algorithm>

namespace opensn
{

struct AAHSweepData
{
  const std::shared_ptr<MeshContinuum>& grid;
  const SpatialDiscretization& discretization;
  const std::vector<UnitCellMatrices>& unit_cell_matrices;
  const std::vector<CellLBSView>& cell_transport_views;
  std::vector<CellOutflowView>& cell_outflow_views;
  const std::vector<double>& source_moments;
  const LBSGroupset& groupset;
  const BlockID2XSMap& xs;
  const unsigned int num_moments;
  const unsigned int max_num_cell_dofs;
  const unsigned int min_num_cell_dofs;
  const bool save_angular_flux;
  const size_t groupset_angle_group_stride;
  const size_t groupset_group_stride;
  std::vector<double>& destination_phi;
  std::vector<double>& destination_psi;
  bool surface_source_active;
  bool include_rhs_time_term;
  DiscreteOrdinatesProblem& problem;
  const std::vector<double>* psi_old; // nullptr for steady sweeps
  unsigned int group_block_size;      // used by fixed-N/AVX path
};

/// Generic sweep kernel (scalar), parameterized by time dependence.
template <bool time_dependent>
inline void
AAH_Sweep_Generic(AAHSweepData& data, AngleSet& angle_set)
{
  const auto& groupset = data.groupset;
  const auto gs_size = groupset.GetNumGroups();
  const auto gs_gi = groupset.first_group;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();

  DenseMatrix<double> Amat(data.max_num_cell_dofs, data.max_num_cell_dofs);
  DenseMatrix<double> Atemp(data.max_num_cell_dofs, data.max_num_cell_dofs);
  std::vector<Vector<double>> b(gs_size, Vector<double>(data.max_num_cell_dofs, 0.0));
  std::vector<double> source(data.max_num_cell_dofs);

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const size_t num_spls = spls.size();

  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    const auto cell_local_id = spls[spls_index];
    auto& cell = data.grid->GetLocalCell(cell_local_id);
    const auto& cell_transport_view = data.cell_transport_views[cell_local_id];
    auto& cell_outflow_view = data.cell_outflow_views[cell_local_id];
    const auto& cell_mapping = data.discretization.GetCellMapping(cell);
    const size_t cell_num_faces = cell.faces.size();
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();

    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    std::vector<double> face_mu_values(cell_num_faces);

    const auto& sigma_t = data.xs.at(cell.block_id)->GetSigmaTotal();

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

    const auto& G = data.unit_cell_matrices[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = data.unit_cell_matrices[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = data.unit_cell_matrices[cell_local_id].intS_shapeI_shapeJ;

    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;
    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();

    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      const auto direction_num = as_angle_indices[as_ss_idx];
      const auto omega = groupset.quadrature->GetOmega(direction_num);
      const auto wt = groupset.quadrature->GetWeight(direction_num);

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
                                          0,
                                          data.surface_source_active);

            if (not psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
          }
        }
      }

      const double* psi_old =
        (time_dependent and data.psi_old)
          ? &(*data.psi_old)[data.discretization.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)]
          : nullptr;
      const double* m2d_row = m2d_op.data() + direction_num * static_cast<size_t>(data.num_moments);
      const double* d2m_row = d2m_op.data() + direction_num * static_cast<size_t>(data.num_moments);

      for (size_t gsg = 0; gsg < gs_size; ++gsg)
      {
        double sigma_tg = sigma_t[gs_gi + gsg];
        if constexpr (time_dependent)
          sigma_tg += tau_gsg[gsg];

        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          double temp_src = 0.0;
          for (unsigned int m = 0; m < data.num_moments; ++m)
          {
            const size_t ir = cell_transport_view.MapDOF(i, m, gs_gi + gsg);
            temp_src += m2d_row[m] * data.source_moments[ir];
          }

          if constexpr (time_dependent)
          {
            const size_t imap =
              i * data.groupset_angle_group_stride + direction_num * data.groupset_group_stride;
            if (data.include_rhs_time_term and psi_old)
              temp_src += tau_gsg[gsg] * psi_old[imap + gsg];
          }

          source[i] = temp_src;
        }

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

      for (unsigned int m = 0; m < data.num_moments; ++m)
      {
        const auto wn_d2m = d2m_row[m];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const size_t ir = cell_transport_view.MapDOF(i, m, gs_gi);
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            data.destination_phi[ir + gsg] += wn_d2m * b[gsg](i);
        }
      }

      if (data.save_angular_flux)
      {
        double* psi_new =
          &data
             .destination_psi[data.discretization.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)];
        double theta = 1.0;
        double inv_theta = 1.0;
        if constexpr (time_dependent)
        {
          theta = data.problem.GetTheta();
          inv_theta = 1.0 / theta;
        }

        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const size_t imap =
            i * data.groupset_angle_group_stride + direction_num * data.groupset_group_stride;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
          {
            const double psi_sol = b[gsg](i);
            if constexpr (time_dependent)
            {
              const double psi_old_val = psi_old ? psi_old[imap + gsg] : 0.0;
              psi_new[imap + gsg] = inv_theta * (psi_sol + (theta - 1.0) * psi_old_val);
            }
            else
            {
              psi_new[imap + gsg] = psi_sol;
            }
          }
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
        const auto& IntF_shapeI = data.unit_cell_matrices[cell_local_id].intS_shapeI[f];

        if (not is_boundary_face and not is_local_face)
          ++deploc_face_counter;

        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);

          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_outflow_view.Add(
              f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg](i) * IntF_shapeI(i));

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

/**
 * Per-cell dispatch sweep kernel (steady-state only). Unlike Sweep_FixedN<N>, which specializes
 * an entire mesh sweep at compile time and requires every cell to share the same node count N in
 * [2, 8], this kernel picks a solve strategy independently for each cell: SimdBatchSolve<N>
 * (compile-time N) when this cell's own node count is in [2, 8], else the runtime-N
 * SimdBatchSolveDynamic, with a scalar GaussEliminateDynamic tail for any group-block
 * remainder that doesn't fill a whole SIMD lane group. This replaces AAH_Sweep_Generic<false> as
 * the fallback used whenever Sweep_FixedN<N>'s whole-mesh restriction doesn't apply -- meshes with
 * mixed cell topologies, or a uniform N outside [2, 8] -- turning what used to be an
 * all-or-nothing option (one stray cell forces the entire mesh to the fully scalar generic kernel)
 * into a per-cell solve decision.
 *
 * FLUDS slot indexing invariant: deploc_face_counter/preloc_face_counter must be incremented in
 * exactly the same per-cell, per-face order as AAH_Sweep_Generic/AAH_Sweep_FixedN, since
 * AAH_FLUDS::SlotDynamics allocated FLUDS slots assuming that exact traversal order.
 */
inline void
AAH_Sweep_Unified(AAHSweepData& data, AngleSet& angle_set)
{
  const auto& groupset = data.groupset;
  const auto gs_size = groupset.GetNumGroups();
  const auto gs_gi = groupset.first_group;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();

  const size_t max_n = data.max_num_cell_dofs;
  std::vector<double> Amat(max_n * max_n, 0.0);
  std::vector<double> mass_matrix(max_n * max_n, 0.0);
  std::vector<double> b(static_cast<size_t>(gs_size) * max_n, 0.0);
  std::vector<double> sigma_block;
  sigma_block.reserve(data.group_block_size);
  std::vector<double> A_scratch(max_n * max_n * Simd::size, 0.0);
  std::vector<double> rhs_scratch(max_n * Simd::size, 0.0);
  std::vector<double> face_mu_values;
  std::vector<std::vector<size_t>> moment_dof_map(data.num_moments);
  for (auto& row : moment_dof_map)
    row.resize(max_n);

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const size_t num_spls = spls.size();

  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    const auto cell_local_id = spls[spls_index];
    auto& cell = data.grid->GetLocalCell(cell_local_id);
    const auto& cell_transport_view = data.cell_transport_views[cell_local_id];
    auto& cell_outflow_view = data.cell_outflow_views[cell_local_id];
    const auto& cell_mapping = data.discretization.GetCellMapping(cell);
    const size_t cell_num_faces = cell.faces.size();
    const size_t N = cell_mapping.GetNumNodes();
    const auto Idx = [N](size_t i, size_t j) { return i * N + j; };

    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    // Reused across cells via resize() rather than a fresh heap allocation every cell.
    face_mu_values.resize(cell_num_faces);

    const auto& sigma_t = cell_transport_view.GetXS().GetSigmaTotal();

    const auto& G = data.unit_cell_matrices[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = data.unit_cell_matrices[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = data.unit_cell_matrices[cell_local_id].intS_shapeI_shapeJ;

    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < N; ++j)
        mass_matrix[Idx(i, j)] = M(i, j);

    for (unsigned int m = 0; m < data.num_moments; ++m)
      for (size_t i = 0; i < N; ++i)
        moment_dof_map[m][i] = cell_transport_view.MapDOF(i, m, gs_gi);

    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;
    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();

    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      const auto direction_num = as_angle_indices[as_ss_idx];
      const auto omega = groupset.quadrature->GetOmega(direction_num);
      const auto wt = groupset.quadrature->GetWeight(direction_num);

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      std::fill_n(b.begin(), static_cast<size_t>(gs_size) * N, 0.0);

      for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
          Amat[Idx(i, j)] = omega.Dot(G(i, j));

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
            Amat[Idx(i, j)] += mu_Nij;

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
                                          data.surface_source_active);

            if (not psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg * N + i] += psi[gsg] * mu_Nij;
          }
        }
      }

      const auto dir_moment_offset =
        static_cast<size_t>(direction_num) * static_cast<size_t>(data.num_moments);
      const double* m2d_row = m2d_op.data() + dir_moment_offset;
      const double* d2m_row = d2m_op.data() + dir_moment_offset;

      for (unsigned int g0 = 0; g0 < gs_size; g0 += data.group_block_size)
      {
        const auto g1 = std::min(g0 + data.group_block_size, gs_size);
        const auto block_len = g1 - g0;
        sigma_block.resize(block_len);

        for (unsigned int gsg = g0; gsg < g1; ++gsg)
        {
          const size_t rel = gsg - g0;
          sigma_block[rel] = sigma_t[gs_gi + gsg];

          double* __restrict bg = &b[static_cast<size_t>(gsg) * N];
          for (unsigned int m = 0; m < data.num_moments; ++m)
          {
            const double w = m2d_row[m];
            for (size_t i = 0; i < N; ++i)
            {
              double value = 0.0;
              for (size_t j = 0; j < N; ++j)
                value +=
                  mass_matrix[Idx(i, j)] * (w * data.source_moments[moment_dof_map[m][j] + gsg]);
              bg[i] += value;
            }
          }
        }

        size_t k = 0;

        double* b_block = &b[static_cast<size_t>(g0) * N];
        detail::DispatchSimdGroupBlock(N,
                                       block_len,
                                       k,
                                       Amat.data(),
                                       mass_matrix.data(),
                                       sigma_block.data(),
                                       b_block,
                                       A_scratch.data(),
                                       rhs_scratch.data());

        for (; k < block_len; ++k)
        {
          const size_t gsg = g0 + k;
          detail::GaussEliminateDynamic(static_cast<int>(N),
                                        Amat.data(),
                                        mass_matrix.data(),
                                        sigma_block[k],
                                        &b[gsg * N],
                                        A_scratch.data());
        }

        for (size_t gsg = g0; gsg < g1; ++gsg)
        {
          const double* __restrict bg = &b[gsg * N];
          for (unsigned int m = 0; m < data.num_moments; ++m)
          {
            const double w = d2m_row[m];
            for (size_t i = 0; i < N; ++i)
              data.destination_phi[moment_dof_map[m][i] + gsg] += w * bg[i];
          }
        }
      }

      if (data.save_angular_flux)
      {
        double* psi_new =
          &data
             .destination_psi[data.discretization.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)];
        for (size_t i = 0; i < N; ++i)
        {
          const size_t imap =
            i * data.groupset_angle_group_stride + direction_num * data.groupset_group_stride;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            psi_new[imap + gsg] = b[gsg * N + i];
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
        const auto& IntF_shapeI = data.unit_cell_matrices[cell_local_id].intS_shapeI[f];

        if (not is_boundary_face and not is_local_face)
          ++deploc_face_counter;

        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);

          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_outflow_view.Add(
              f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg * N + i] * IntF_shapeI(i));

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
              psi[gsg] = b[gsg * N + i];
          }
        }
      }
    }
  }
}

// Fixed-N kernel (scalar/AVX) is specialized in aah_avx_sweep_chunk.cc.
template <unsigned int NumNodes, bool time_dependent>
void AAH_Sweep_FixedN(AAHSweepData& data, AngleSet& angle_set);

} // namespace opensn
