// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"

#if defined(__AVX512F__)
#include <immintrin.h>

namespace opensn
{

// Solve 8 independent 4x4 systems in SIMD: (Amat + sigma*M) * x = b
// Inputs:
//   Am[16]    : base coefficients of Amat (row-major: 00,01,...,33)
//   Mm[16]    : base coefficients of M (same layout)
//   sigma[8]  : sigma_tg for 8 consecutive groups
//   b_base    : pointer into AoS b array at the first group in the batch
//               layout is [g0*4+0, g0*4+1, g0*4+2, g0*4+3, g1*4+0, ...]
// Effects:
//   Overwrites the corresponding 8 b-vectors in place with the solutions.
static inline void
BatchSolveAVX512(const double* Am, const double* Mm, const double* sigma, double* b_base)
{
  // Build static byte offsets (relative to b_base) for gathers/scatters.
  // AoS layout: for lane k in [0..7], element i in [0..3] is at offset (4*k + i)
  const __m512i idx0 =
    _mm512_setr_epi64(0 * 8, 4 * 8, 8 * 8, 12 * 8, 16 * 8, 20 * 8, 24 * 8, 28 * 8);
  const __m512i idx1 =
    _mm512_setr_epi64(1 * 8, 5 * 8, 9 * 8, 13 * 8, 17 * 8, 21 * 8, 25 * 8, 29 * 8);
  const __m512i idx2 =
    _mm512_setr_epi64(2 * 8, 6 * 8, 10 * 8, 14 * 8, 18 * 8, 22 * 8, 26 * 8, 30 * 8);
  const __m512i idx3 =
    _mm512_setr_epi64(3 * 8, 7 * 8, 11 * 8, 15 * 8, 19 * 8, 23 * 8, 27 * 8, 31 * 8);

  // Load RHS batch
  char* bbytes = reinterpret_cast<char*>(b_base);
  __m512d b0 = _mm512_i64gather_pd(idx0, bbytes, 1);
  __m512d b1 = _mm512_i64gather_pd(idx1, bbytes, 1);
  __m512d b2 = _mm512_i64gather_pd(idx2, bbytes, 1);
  __m512d b3 = _mm512_i64gather_pd(idx3, bbytes, 1);

  // Load sigma
  __m512d sv = _mm512_loadu_pd(sigma);

  // Broadcast Am and M scalars
  auto bcast = [](double x) { return _mm512_set1_pd(x); };

  const __m512d Am00 = bcast(Am[0]), Am01 = bcast(Am[1]), Am02 = bcast(Am[2]), Am03 = bcast(Am[3]);
  const __m512d Am10 = bcast(Am[4]), Am11 = bcast(Am[5]), Am12 = bcast(Am[6]), Am13 = bcast(Am[7]);
  const __m512d Am20 = bcast(Am[8]), Am21 = bcast(Am[9]), Am22 = bcast(Am[10]),
                Am23 = bcast(Am[11]);
  const __m512d Am30 = bcast(Am[12]), Am31 = bcast(Am[13]), Am32 = bcast(Am[14]),
                Am33 = bcast(Am[15]);

  const __m512d M00 = bcast(Mm[0]), M01 = bcast(Mm[1]), M02 = bcast(Mm[2]), M03 = bcast(Mm[3]);
  const __m512d M10 = bcast(Mm[4]), M11 = bcast(Mm[5]), M12 = bcast(Mm[6]), M13 = bcast(Mm[7]);
  const __m512d M20 = bcast(Mm[8]), M21 = bcast(Mm[9]), M22 = bcast(Mm[10]), M23 = bcast(Mm[11]);
  const __m512d M30 = bcast(Mm[12]), M31 = bcast(Mm[13]), M32 = bcast(Mm[14]), M33 = bcast(Mm[15]);

  // A = Am + sigma*M (elementwise)
  auto madd = [&sv](const __m512d& a, const __m512d& m) { return _mm512_fmadd_pd(sv, m, a); };

  __m512d A00 = madd(Am00, M00), A01 = madd(Am01, M01), A02 = madd(Am02, M02),
          A03 = madd(Am03, M03);
  __m512d A10 = madd(Am10, M10), A11 = madd(Am11, M11), A12 = madd(Am12, M12),
          A13 = madd(Am13, M13);
  __m512d A20 = madd(Am20, M20), A21 = madd(Am21, M21), A22 = madd(Am22, M22),
          A23 = madd(Am23, M23);
  __m512d A30 = madd(Am30, M30), A31 = madd(Am31, M31), A32 = madd(Am32, M32),
          A33 = madd(Am33, M33);

  // Forward elimination (SIMD lanes are independent systems)
  const __m512d invA00 = _mm512_div_pd(_mm512_set1_pd(1.0), A00);

  const __m512d v10 = _mm512_mul_pd(A10, invA00);
  b1 = _mm512_fnmadd_pd(v10, b0, b1);
  A11 = _mm512_fnmadd_pd(v10, A01, A11);
  A12 = _mm512_fnmadd_pd(v10, A02, A12);
  A13 = _mm512_fnmadd_pd(v10, A03, A13);

  const __m512d v20 = _mm512_mul_pd(A20, invA00);
  b2 = _mm512_fnmadd_pd(v20, b0, b2);
  A21 = _mm512_fnmadd_pd(v20, A01, A21);
  A22 = _mm512_fnmadd_pd(v20, A02, A22);
  A23 = _mm512_fnmadd_pd(v20, A03, A23);

  const __m512d v30 = _mm512_mul_pd(A30, invA00);
  b3 = _mm512_fnmadd_pd(v30, b0, b3);
  A31 = _mm512_fnmadd_pd(v30, A01, A31);
  A32 = _mm512_fnmadd_pd(v30, A02, A32);
  A33 = _mm512_fnmadd_pd(v30, A03, A33);

  const __m512d invA11 = _mm512_div_pd(_mm512_set1_pd(1.0), A11);

  const __m512d v21 = _mm512_mul_pd(A21, invA11);
  b2 = _mm512_fnmadd_pd(v21, b1, b2);
  A22 = _mm512_fnmadd_pd(v21, A12, A22);
  A23 = _mm512_fnmadd_pd(v21, A13, A23);

  const __m512d v31 = _mm512_mul_pd(A31, invA11);
  b3 = _mm512_fnmadd_pd(v31, b1, b3);
  A32 = _mm512_fnmadd_pd(v31, A12, A32);
  A33 = _mm512_fnmadd_pd(v31, A13, A33);

  const __m512d invA22 = _mm512_div_pd(_mm512_set1_pd(1.0), A22);

  const __m512d v32 = _mm512_mul_pd(A32, invA22);
  b3 = _mm512_fnmadd_pd(v32, b2, b3);
  A33 = _mm512_fnmadd_pd(v32, A23, A33);

  // Back substitution
  b3 = _mm512_div_pd(b3, A33);
  b2 = _mm512_mul_pd(_mm512_sub_pd(b2, _mm512_mul_pd(A23, b3)), invA22);
  b1 = _mm512_mul_pd(
    _mm512_sub_pd(_mm512_sub_pd(b1, _mm512_mul_pd(A12, b2)), _mm512_mul_pd(A13, b3)), invA11);
  b0 = _mm512_mul_pd(
    _mm512_sub_pd(_mm512_sub_pd(_mm512_sub_pd(b0, _mm512_mul_pd(A01, b1)), _mm512_mul_pd(A02, b2)),
                  _mm512_mul_pd(A03, b3)),
    invA00);

  // Scatter solutions back to AoS
  _mm512_i64scatter_pd(bbytes, idx0, b0, 1);
  _mm512_i64scatter_pd(bbytes, idx1, b1, 1);
  _mm512_i64scatter_pd(bbytes, idx2, b2, 1);
  _mm512_i64scatter_pd(bbytes, idx3, b3, 1);
}
#endif

void
AAHSweepChunk::CPUSweep_N4(AngleSet& angle_set)
{
  auto gs_size = groupset_.groups.size();
  auto gs_gi = groupset_.groups.front().id;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();

  DenseMatrix<double> Amat(4, 4);
  DenseMatrix<double> Atemp(4, 4);
  std::vector<double> b(groupset_.groups.size() * 4, 0.0);
  // std::vector<Vector<double>> b(groupset_.groups.size(), Vector<double>(4, 0.0));
  std::vector<double> source(4);

  // Loop over each cell
  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const size_t num_spls = spls.size();
  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    auto cell_local_id = spls[spls_index];
    auto& cell = grid_->local_cells[cell_local_id];
    auto& cell_transport_view = cell_transport_views_[cell_local_id];
    const auto& cell_mapping = discretization_.GetCellMapping(cell);

    // Cell face data
    std::vector<double> face_mu_values(4);
    const auto cell_num_faces = cell.faces.size();
    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;

    // Material properties
    const auto& rho = densities_[cell.local_id];
    const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();

    // Cell matrices
    const auto& G = unit_cell_matrices_[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = unit_cell_matrices_[cell_local_id].intS_shapeI_shapeJ;

    // Loop over angles in angleset (as = angleset, ss = subset)
    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      auto direction_num = as_angle_indices[as_ss_idx];
      auto omega = groupset_.quadrature->omegas[direction_num];
      auto wt = groupset_.quadrature->weights[direction_num];

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      // Reset right-hand side
      for (size_t gsg = 0; gsg < gs_size; ++gsg)
        for (int i = 0; i < 4; ++i)
          b[gsg * 4 + i] = 0.0;

      // Initialize A matrix
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          Amat(i, j) = omega.Dot(G(i, j));

      // Update face orientations
      for (int f = 0; f < 4; ++f)
        face_mu_values[f] = omega.Dot(cell.faces[f].normal);

      // Surface integrals
      int in_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::INCOMING)
          continue;

        auto& cell_face = cell.faces[f];
        const auto& Ms_f = M_surf[f];
        const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary_face = !cell_face.has_neighbor;
        const double mu_f = -face_mu_values[f];

        if (is_local_face)
          ++in_face_counter;
        else if (!is_boundary_face)
          ++preloc_face_counter;

        for (size_t fj = 0; fj < num_face_nodes; ++fj)
        {
          const int j = cell_mapping.MapFaceNode(f, fj);
          const double* psi = nullptr;

          if (is_local_face)
            psi = fluds.UpwindPsi(spls_index, in_face_counter, fj, 0, as_ss_idx);
          else if (!is_boundary_face)
            psi = fluds.NLUpwindPsi(preloc_face_counter, fj, 0, as_ss_idx);
          else
            psi = angle_set.PsiBoundary(cell_face.neighbor_id,
                                        direction_num,
                                        cell_local_id,
                                        f,
                                        fj,
                                        gs_gi,
                                        IsSurfaceSourceActive());

          for (size_t fi = 0; fi < num_face_nodes; ++fi)
          {
            const int i = cell_mapping.MapFaceNode(f, fi);
            const double mu_Nij = mu_f * Ms_f(i, j);

            Amat(i, j) += mu_Nij;

            if (!psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg * 4 + i] += psi[gsg] * mu_Nij;
          }
        }
      }

      const double M00 = M(0, 0), M01 = M(0, 1), M02 = M(0, 2), M03 = M(0, 3);
      const double M10 = M(1, 0), M11 = M(1, 1), M12 = M(1, 2), M13 = M(1, 3);
      const double M20 = M(2, 0), M21 = M(2, 1), M22 = M(2, 2), M23 = M(2, 3);
      const double M30 = M(3, 0), M31 = M(3, 1), M32 = M(3, 2), M33 = M(3, 3);

      const double Mm[16] = {
        M00, M01, M02, M03, M10, M11, M12, M13, M20, M21, M22, M23, M30, M31, M32, M33};

      const double Am00 = Amat(0, 0), Am01 = Amat(0, 1), Am02 = Amat(0, 2), Am03 = Amat(0, 3);
      const double Am10 = Amat(1, 0), Am11 = Amat(1, 1), Am12 = Amat(1, 2), Am13 = Amat(1, 3);
      const double Am20 = Amat(2, 0), Am21 = Amat(2, 1), Am22 = Amat(2, 2), Am23 = Amat(2, 3);
      const double Am30 = Amat(3, 0), Am31 = Amat(3, 1), Am32 = Amat(3, 2), Am33 = Amat(3, 3);

      const double Am[16] = {Am00,
                             Am01,
                             Am02,
                             Am03,
                             Am10,
                             Am11,
                             Am12,
                             Am13,
                             Am20,
                             Am21,
                             Am22,
                             Am23,
                             Am30,
                             Am31,
                             Am32,
                             Am33};

      const double* __restrict m2d_row = m2d_op[direction_num].data();
      const double* __restrict d2m_row = d2m_op[direction_num].data();

      // Process groups in blocks (keep your group_block_size_, but try to make it a multiple of 8)
      for (size_t g0 = 0; g0 < gs_size; g0 += group_block_size_)
      {
        const size_t g1 = std::min(g0 + group_block_size_, gs_size);

        // Phase 1: assemble RHS b for this block (and cache sigma_tg)
        // We also drop 4x MapDOF calls per group per moment by using base+gsg.
        std::vector<double> sigma_blk(g1 - g0);

        for (size_t gsg = g0; gsg < g1; ++gsg)
        {
          const size_t rel = gsg - g0;
          const double sigma_tg = rho * sigma_t[gs_gi + gsg];
          sigma_blk[rel] = sigma_tg;

          // Add source moments: b_i += sum_m M(i, :) * (w * q_moms at nodes 0..3)
          for (int m = 0; m < num_moments_; ++m)
          {
            const double w = m2d_row[m];

            // Base indices for this (cell, m, group-set start)
            const size_t ir0_base = cell_transport_view.MapDOF(0, m, gs_gi);
            const size_t ir1_base = cell_transport_view.MapDOF(1, m, gs_gi);
            const size_t ir2_base = cell_transport_view.MapDOF(2, m, gs_gi);
            const size_t ir3_base = cell_transport_view.MapDOF(3, m, gs_gi);

            const double s0 = w * source_moments_[ir0_base + gsg];
            const double s1 = w * source_moments_[ir1_base + gsg];
            const double s2 = w * source_moments_[ir2_base + gsg];
            const double s3 = w * source_moments_[ir3_base + gsg];

            // Accumulate with scalar FMAs; the compiler will fuse
            double* __restrict bg = &b[gsg * 4];
            bg[0] += M00 * s0 + M01 * s1 + M02 * s2 + M03 * s3;
            bg[1] += M10 * s0 + M11 * s1 + M12 * s2 + M13 * s3;
            bg[2] += M20 * s0 + M21 * s1 + M22 * s2 + M23 * s3;
            bg[3] += M30 * s0 + M31 * s1 + M32 * s2 + M33 * s3;
          }
        }

        // Phase 2: solve the 4x4 systems for this block
        size_t blk_len = g1 - g0;
        size_t k = 0;

#if defined(__AVX512F__)
        // Vectorized batches of 8 groups
        for (; k + 8 <= blk_len; k += 8)
          BatchSolveAVX512(Am, Mm, &sigma_blk[k], &b[(g0 + k) * 4]);
#endif

        // Remainder (or entire block if no AVX-512)
        for (; k < blk_len; ++k)
        {
          const size_t gsg = g0 + k;
          const double sigma_tg = sigma_blk[k];

          // Build Atemp = Amat + sigma_tg*M
          double A00 = Am00 + sigma_tg * M00, A01 = Am01 + sigma_tg * M01;
          double A02 = Am02 + sigma_tg * M02, A03 = Am03 + sigma_tg * M03;

          double A10 = Am10 + sigma_tg * M10, A11 = Am11 + sigma_tg * M11;
          double A12 = Am12 + sigma_tg * M12, A13 = Am13 + sigma_tg * M13;

          double A20 = Am20 + sigma_tg * M20, A21 = Am21 + sigma_tg * M21;
          double A22 = Am22 + sigma_tg * M22, A23 = Am23 + sigma_tg * M23;

          double A30 = Am30 + sigma_tg * M30, A31 = Am31 + sigma_tg * M31;
          double A32 = Am32 + sigma_tg * M32, A33 = Am33 + sigma_tg * M33;

          double* __restrict bg = &b[gsg * 4];

          // Forward elimination (scalar fallback)
          const double invA00 = 1.0 / A00;

          const double v10 = A10 * invA00;
          bg[1] -= v10 * bg[0];
          A11 -= v10 * A01;
          A12 -= v10 * A02;
          A13 -= v10 * A03;

          const double v20 = A20 * invA00;
          bg[2] -= v20 * bg[0];
          A21 -= v20 * A01;
          A22 -= v20 * A02;
          A23 -= v20 * A03;

          const double v30 = A30 * invA00;
          bg[3] -= v30 * bg[0];
          A31 -= v30 * A01;
          A32 -= v30 * A02;
          A33 -= v30 * A03;

          const double invA11 = 1.0 / A11;

          const double v21 = A21 * invA11;
          bg[2] -= v21 * bg[1];
          A22 -= v21 * A12;
          A23 -= v21 * A13;

          const double v31 = A31 * invA11;
          bg[3] -= v31 * bg[1];
          A32 -= v31 * A12;
          A33 -= v31 * A13;

          const double invA22 = 1.0 / A22;

          const double v32 = A32 * invA22;
          bg[3] -= v32 * bg[2];
          A33 -= v32 * A23;

          // Back substitution
          bg[3] = bg[3] / A33;
          bg[2] = (bg[2] - A23 * bg[3]) * invA22;
          bg[1] = (bg[1] - A12 * bg[2] - A13 * bg[3]) * invA11;
          bg[0] = (bg[0] - A01 * bg[1] - A02 * bg[2] - A03 * bg[3]) * invA00;
        }

        // Phase 3: moment accumulation into phi (unchanged math, but hoisted row ptrs)
        for (size_t gsg = g0; gsg < g1; ++gsg)
        {
          const double* __restrict bg = &b[gsg * 4];

          for (int m = 0; m < num_moments_; ++m)
          {
            const double w = d2m_row[m];

            const size_t ir0_base = cell_transport_view.MapDOF(0, m, gs_gi);
            const size_t ir1_base = cell_transport_view.MapDOF(1, m, gs_gi);
            const size_t ir2_base = cell_transport_view.MapDOF(2, m, gs_gi);
            const size_t ir3_base = cell_transport_view.MapDOF(3, m, gs_gi);

            destination_phi_[ir0_base + gsg] += w * bg[0];
            destination_phi_[ir1_base + gsg] += w * bg[1];
            destination_phi_[ir2_base + gsg] += w * bg[2];
            destination_phi_[ir3_base + gsg] += w * bg[3];
          }
        }
      } // blocks

      // Save angular flux during sweep
      if (save_angular_flux_)
      {
        double* cell_psi_data =
          &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];

        for (int i = 0; i < 4; ++i)
        {
          const size_t imap =
            i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_psi_data[imap + gsg] = b[gsg * 4 + i];
        }
      }

      // For outoing, non-boundary faces, copy angular flux to fluds and
      // accumulate outflow
      int out_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::OUTGOING)
          continue;

        const auto& face = cell.faces[f];
        const auto& IntF_shapeI = unit_cell_matrices_[cell_local_id].intS_shapeI[f];
        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary = !face.has_neighbor;
        const bool is_reflecting =
          is_boundary && angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting();
        const double mu_wt_f = wt * face_mu_values[f];

        ++out_face_counter;

        if (!is_boundary && !is_local_face)
          ++deploc_face_counter;

        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);

          if (is_boundary)
          {
            const double flux_i = mu_wt_f * IntF_shapeI(i);

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              cell_transport_view.AddOutflow(f, gs_gi + gsg, flux_i * b[gsg * 4 + i]);
          }

          double* psi = nullptr;
          if (is_local_face)
            psi = fluds.OutgoingPsi(spls_index, out_face_counter, fi, as_ss_idx);
          else if (!is_boundary)
            psi = fluds.NLOutgoingPsi(deploc_face_counter, fi, as_ss_idx);
          else if (is_reflecting)
            psi = angle_set.PsiReflected(face.neighbor_id, direction_num, cell_local_id, f, fi);
          else
            continue; // non-reflecting boundary has no psi target

          // Write psi only for interior or reflecting boundary faces
          if (!is_boundary || is_reflecting)
          {
            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              psi[gsg] = b[gsg * 4 + i];
          }
        }
      }
    } // for angleset/subset
  }   // for cell
}

} // namespace opensn
