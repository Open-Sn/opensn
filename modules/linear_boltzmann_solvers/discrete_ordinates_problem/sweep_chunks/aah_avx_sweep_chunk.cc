// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log_exceptions.h"
#include "caliper/cali.h"
#include <algorithm>
#include <array>
#include <vector>

#if __AVX512F__ || __AVX2__
#include <immintrin.h>
#endif

#if __clang__ || __INTEL_COMPILER
#define PRAGMA_UNROLL _Pragma("unroll")
#elif __GNUC__
#define PRAGMA_UNROLL _Pragma("GCC unroll 8")
#else
#define PRAGMA_UNROLL
#endif

namespace opensn
{

namespace detail
{

#if __AVX512F__
struct AVX512Ops
{
  using avx_vec = __m512d;
  using avx_index = __m512i;

  static inline avx_vec LoadSigma(const double* sigma) { return _mm512_loadu_pd(sigma); }
  static inline avx_vec Set1(double x) { return _mm512_set1_pd(x); }
  static inline avx_vec Add(const avx_vec& a, const avx_vec& b) { return _mm512_add_pd(a, b); }
  static inline avx_vec Sub(const avx_vec& a, const avx_vec& b) { return _mm512_sub_pd(a, b); }
  static inline avx_vec Mul(const avx_vec& a, const avx_vec& b) { return _mm512_mul_pd(a, b); }
  static inline avx_vec Div(const avx_vec& a, const avx_vec& b) { return _mm512_div_pd(a, b); }
  static inline avx_vec Fmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    // a + b * c
    return _mm512_fmadd_pd(b, c, a);
  }
  static inline avx_vec Fnmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    // c - a * b
    return _mm512_fnmadd_pd(a, b, c);
  }
  static inline avx_vec Reciprocal(const avx_vec& v) { return Div(Set1(1.0), v); }
  static inline avx_vec Gather(const avx_index& idx, const double* base)
  {
    return _mm512_i64gather_pd(idx, base, sizeof(double));
  }
  static inline void Scatter(const avx_index& idx, double* base, const avx_vec& value)
  {
    _mm512_i64scatter_pd(base, idx, value, sizeof(double));
  }
};
#elif __AVX2__
struct AVX2Ops
{
  using avx_vec = __m256d;
  using avx_index = __m128i;

  static inline avx_vec LoadSigma(const double* sigma) { return _mm256_loadu_pd(sigma); }
  static inline avx_vec Set1(double x) { return _mm256_set1_pd(x); }
  static inline avx_vec Add(const avx_vec& a, const avx_vec& b) { return _mm256_add_pd(a, b); }
  static inline avx_vec Sub(const avx_vec& a, const avx_vec& b) { return _mm256_sub_pd(a, b); }
  static inline avx_vec Mul(const avx_vec& a, const avx_vec& b) { return _mm256_mul_pd(a, b); }
  static inline avx_vec Div(const avx_vec& a, const avx_vec& b) { return _mm256_div_pd(a, b); }

#if __FMA__
  static inline avx_vec Fmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    return _mm256_fmadd_pd(b, c, a);
  }
  static inline avx_vec Fnmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    return _mm256_fnmadd_pd(a, b, c);
  }
#else
  static inline Vec Fmadd(const Vec& a, const Vec& b, const Vec& c) { return Add(a, Mul(b, c)); }
  static inline Vec Fnmadd(const Vec& a, const Vec& b, const Vec& c) { return Sub(c, Mul(a, b)); }
#endif

  static inline avx_vec Reciprocal(const avx_vec& v) { return Div(Set1(1.0), v); }
  static inline avx_vec Gather(const avx_index& idx, const double* base)
  {
    return _mm256_i32gather_pd(base, idx, sizeof(double));
  }
  static inline void Scatter(const avx_index& idx, double* base, const avx_vec& value)
  {
    alignas(32) double buffer[simd_width];
    _mm256_store_pd(buffer, value);
    alignas(16) int offsets[simd_width];
    _mm_store_si128(reinterpret_cast<__m128i*>(offsets), idx);
    for (int lane = 0; lane < simd_width; ++lane)
      base[offsets[lane]] = buffer[lane];
  }
};
#endif

template <class Ops, int N>
struct GatherIndexBuilder
{
  static typename Ops::avx_index Build(int /*unused*/)
  {
    static_assert(sizeof(Ops) == 0, "SIMD gather index helper not implemented for this Ops type.");
    return typename Ops::avx_index{};
  }
};

#if __AVX512F__
template <int N>
struct GatherIndexBuilder<AVX512Ops, N>
{
  static AVX512Ops::avx_index Build(int row)
  {
    long long vals[simd_width];
    for (int lane = 0; lane < simd_width; ++lane)
      vals[lane] = static_cast<long long>(lane * N + row);
    return _mm512_setr_epi64(
      vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7]);
  }
};
#elif __AVX2__
template <int N>
struct GatherIndexBuilder<AVX2Ops, N>
{
  static AVX2Ops::avx_index Build(int row)
  {
    int vals[simd_width];
    for (int lane = 0; lane < simd_width; ++lane)
      vals[lane] = lane * N + row;
    return _mm_setr_epi32(vals[0], vals[1], vals[2], vals[3]);
  }
};
#endif

template <class Ops, int N>
inline typename Ops::avx_index
MakeGatherIndex(int row)
{
  return GatherIndexBuilder<Ops, N>::Build(row);
}

template <class Ops, int N>
inline void
SimdBatchSolve(const double* Am, const double* Mm, const double* sigma_t, double* __restrict b)
{
  using avx_vec = typename Ops::avx_vec;

  avx_vec rhs[N];
  PRAGMA_UNROLL
  for (int row = 0; row < N; ++row)
    rhs[row] = Ops::Gather(MakeGatherIndex<Ops, N>(row), b);

  const avx_vec sigma = Ops::LoadSigma(sigma_t);
  avx_vec A[N * N];
  PRAGMA_UNROLL
  for (int i = 0; i < N; ++i)
  {
    PRAGMA_UNROLL
    for (int j = 0; j < N; ++j)
    {
      const avx_vec Amij = Ops::Set1(Am[i * N + j]);
      const avx_vec Mmij = Ops::Set1(Mm[i * N + j]);
      A[i * N + j] = Ops::Fmadd(Amij, sigma, Mmij);
    }
  }

  auto entry = [&](int i, int j) -> avx_vec& { return A[i * N + j]; };
  PRAGMA_UNROLL
  for (int pivot = 0; pivot < N; ++pivot)
  {
    const avx_vec inv = Ops::Reciprocal(entry(pivot, pivot));
    PRAGMA_UNROLL
    for (int row = pivot + 1; row < N; ++row)
    {
      const avx_vec factor = Ops::Mul(entry(row, pivot), inv);
      rhs[row] = Ops::Fnmadd(factor, rhs[pivot], rhs[row]);
      PRAGMA_UNROLL
      for (int col = pivot + 1; col < N; ++col)
        entry(row, col) = Ops::Fnmadd(factor, entry(pivot, col), entry(row, col));
    }
  }

  PRAGMA_UNROLL
  for (int pivot = N - 1; pivot >= 0; --pivot)
  {
    avx_vec rhs_vec = rhs[pivot];
    PRAGMA_UNROLL
    for (int col = pivot + 1; col < N; ++col)
      rhs_vec = Ops::Fnmadd(entry(pivot, col), rhs[col], rhs_vec);
    rhs[pivot] = Ops::Mul(rhs_vec, Ops::Reciprocal(entry(pivot, pivot)));
  }

  PRAGMA_UNROLL
  for (int row = 0; row < N; ++row)
    Ops::Scatter(MakeGatherIndex<Ops, N>(row), b, rhs[row]);
}

} // namespace detail

template <int NumNodes>
void
AAHSweepChunk::CPUSweep_FixedN(AngleSet& angle_set)
{
  static_assert(NumNodes >= 1 and NumNodes <= 8);

  CALI_CXX_MARK_SCOPE("AAHSweepChunk::CPUSweep_FixedN");

  const size_t gs_size = groupset_.groups.size();
  const size_t gs_gi = groupset_.groups.front().id;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();

  std::vector<double> b(gs_size * NumNodes, 0.0);
  std::vector<double> sigma_block;
  sigma_block.reserve(group_block_size_);

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const size_t num_spls = spls.size();
  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    const uint64_t cell_local_id = spls[spls_index];
    auto& cell = grid_->local_cells[cell_local_id];
    auto& cell_transport_view = cell_transport_views_[cell_local_id];
    const auto& cell_mapping = discretization_.GetCellMapping(cell);
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();
    constexpr auto expected_nodes = static_cast<size_t>(NumNodes);
    OpenSnInvalidArgumentIf(
      cell_num_nodes != expected_nodes,
      "AAHSweepChunk::CPUSweep_FixedN invoked for an incompatible cell topology.");

    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    const size_t cell_num_faces = cell.faces.size();
    std::vector<double> face_mu_values(cell_num_faces, 0.0);

    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;

    const double rho = densities_[cell.local_id];
    const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();

    const auto& unit_mats = unit_cell_matrices_[cell_local_id];
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
    std::vector<std::array<size_t, NumNodes>> moment_dof_map(num_moments_);
    for (std::size_t m = 0; m < num_moments_; ++m)
    {
      PRAGMA_UNROLL
      for (int i = 0; i < NumNodes; ++i)
        moment_dof_map[m][i] = cell_transport_view.MapDOF(i, m, gs_gi);
    }

    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      const auto direction_num = as_angle_indices[as_ss_idx];
      const auto omega = groupset_.quadrature->omegas[direction_num];
      const auto wt = groupset_.quadrature->weights[direction_num];

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
                                        static_cast<int>(gs_gi),
                                        IsSurfaceSourceActive());

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

      const double* __restrict m2d_row = m2d_op[direction_num].data();
      const double* __restrict d2m_row = d2m_op[direction_num].data();

      for (size_t g0 = 0; g0 < gs_size; g0 += group_block_size_)
      {
        const size_t g1 = std::min(g0 + group_block_size_, gs_size);
        const size_t block_len = g1 - g0;
        sigma_block.resize(block_len);

        for (size_t gsg = g0; gsg < g1; ++gsg)
        {
          const size_t rel = gsg - g0;
          const double sigma_tg = rho * sigma_t[gs_gi + gsg];
          sigma_block[rel] = sigma_tg;

          double* __restrict bg = &b[gsg * NumNodes];
          for (std::size_t m = 0; m < num_moments_; ++m)
          {
            const double w = m2d_row[m];
            std::array<double, NumNodes> nodal_source{};
            for (int i = 0; i < NumNodes; ++i)
              nodal_source[i] = w * source_moments_[moment_dof_map[m][i] + gsg];

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
          for (std::size_t m = 0; m < num_moments_; ++m)
          {
            const double w = d2m_row[m];
            PRAGMA_UNROLL
            for (int i = 0; i < NumNodes; ++i)
            {
              const size_t dof = cell_transport_view.MapDOF(i, m, gs_gi);
              destination_phi_[dof + gsg] += w * bg[i];
            }
          }
        }
      }

      if (save_angular_flux_)
      {
        double* cell_psi_data =
          &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
        PRAGMA_UNROLL
        for (int i = 0; i < NumNodes; ++i)
        {
          const size_t imap =
            i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_psi_data[imap + gsg] = b[gsg * NumNodes + i];
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

template void AAHSweepChunk::CPUSweep_FixedN<2>(AngleSet&);
template void AAHSweepChunk::CPUSweep_FixedN<3>(AngleSet&);
template void AAHSweepChunk::CPUSweep_FixedN<4>(AngleSet&);
template void AAHSweepChunk::CPUSweep_FixedN<5>(AngleSet&);
template void AAHSweepChunk::CPUSweep_FixedN<6>(AngleSet&);
template void AAHSweepChunk::CPUSweep_FixedN<7>(AngleSet&);
template void AAHSweepChunk::CPUSweep_FixedN<8>(AngleSet&);

} // namespace opensn
