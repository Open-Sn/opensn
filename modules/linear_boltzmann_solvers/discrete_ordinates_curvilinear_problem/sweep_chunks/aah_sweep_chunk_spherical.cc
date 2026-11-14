// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/sweep_chunks/aah_sweep_chunk_spherical.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/discrete_ordinates_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/mesh/mesh/mesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace opensn
{
namespace
{

const CurvilinearProductQuadrature&
RequireSphericalQuadrature(const LBSGroupset& groupset)
{
  const auto* quadrature =
    dynamic_cast<const GLProductQuadrature1DSpherical*>(groupset.quadrature.get());
  if (quadrature == nullptr)
    throw std::invalid_argument(
      "AAHSweepChunkSpherical: 1D spherical geometry requires a GLProductQuadrature1DSpherical.");
  return *quadrature;
}

} // namespace

AAHSweepChunkSpherical::AAHSweepChunkSpherical(DiscreteOrdinatesProblem& problem,
                                               LBSGroupset& groupset)
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
    curvilinear_quadrature_(RequireSphericalQuadrature(groupset)),
    group_block_size_(ComputeGroupBlockSize(groupset.GetNumGroups())),
    secondary_unit_cell_matrices_(dynamic_cast<const DiscreteOrdinatesCurvilinearProblem&>(problem)
                                    .GetSecondaryUnitCellMatrices())
{
  unknown_manager_.AddUnknown(UnknownType::VECTOR_N, groupset_.GetNumGroups());
  psi_edge_.assign(discretization_.GetNumLocalDOFs(unknown_manager_), 0.0);
}

void
AAHSweepChunkSpherical::Sweep(AngleSet& angle_set)
{
  // PWLD radial slab cells have two nodes. Retain the runtime-size fallback.
  if (min_num_cell_dofs_ == 2 and max_num_cell_dofs_ == 2)
    Sweep_FixedN<2>(angle_set);
  else
    Sweep_Generic(angle_set);
}

void
AAHSweepChunkSpherical::Sweep_Generic(AngleSet& angle_set)
{
  const auto gs_size = groupset_.GetNumGroups();
  const auto gs_gi = groupset_.first_group;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();
  const auto& fac_diamond_difference = curvilinear_quadrature_.GetDiamondDifferenceFactor();
  const auto& fac_streaming_operator = curvilinear_quadrature_.GetStreamingOperatorFactor();

  DenseMatrix<double> Amat(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Atemp(max_num_cell_dofs_, max_num_cell_dofs_);
  std::vector<Vector<double>> b(gs_size, Vector<double>(max_num_cell_dofs_));
  std::vector<double> source(max_num_cell_dofs_);

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const auto& as_angle_indices = angle_set.GetAngleIndices();
  for (size_t spls_index = 0; spls_index < spls.size(); ++spls_index)
  {
    const auto cell_local_id = spls[spls_index];
    const auto& cell = grid_->GetLocalCell(cell_local_id);
    const auto& cell_mapping = discretization_.GetCellMapping(cell);
    const auto& cell_transport_view = cell_transport_views_[cell_local_id];
    auto& cell_outflow_view = cell_outflow_views_[cell_local_id];
    const auto cell_num_faces = cell.faces.size();
    const auto cell_num_nodes = cell_mapping.GetNumNodes();
    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    std::vector<double> face_mu_values(cell_num_faces);

    const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();
    const auto& G = unit_cell_matrices_[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = unit_cell_matrices_[cell_local_id].intS_shapeI_shapeJ;
    const auto& M_r = secondary_unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;

    const auto ni_deploc_face_counter = deploc_face_counter;
    const auto ni_preloc_face_counter = preloc_face_counter;
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      const auto direction_num = as_angle_indices[as_ss_idx];
      const auto& omega = groupset_.quadrature->GetOmega(direction_num);
      const auto wt = groupset_.quadrature->GetWeight(direction_num);
      const double tau = fac_diamond_difference[direction_num];
      const double c_ang = fac_streaming_operator[direction_num];
      // The zero-weight final direction mu = +1 ends the recursion. Its flux is the angular-edge
      // flux psi_{N+1/2}; solving for it would be singular in a void cell at the center, whose
      // inflow face has zero area. It contributes nothing to the flux moments and is used only as
      // reflected inflow for the starting direction.
      const bool final_direction = wt == 0.0 and omega.z > 0.0;

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      // Angular-redistribution term: c_ang (psi_n - psi_{n-1/2}) / r
      for (size_t gsg = 0; gsg < gs_size; ++gsg)
        b[gsg] = Vector<double>(cell_num_nodes, 0.0);
      for (size_t i = 0; i < cell_num_nodes; ++i)
        for (size_t j = 0; j < cell_num_nodes; ++j)
        {
          const auto jr = discretization_.MapDOFLocal(cell, j, unknown_manager_, 0, 0);
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            b[gsg](i) += c_ang * M_r(i, j) * psi_edge_[jr + gsg];
        }

      for (size_t i = 0; i < cell_num_nodes; ++i)
        for (size_t j = 0; j < cell_num_nodes; ++j)
          Amat(i, j) = omega.Dot(G(i, j)) + c_ang * M_r(i, j);

      for (size_t f = 0; f < cell_num_faces; ++f)
        face_mu_values[f] = omega.Dot(cell.faces[f].normal);

      // Incoming surface terms
      int in_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::INCOMING)
          continue;

        const auto& cell_face = cell.faces[f];
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary_face = not cell_face.has_neighbor;
        // The center of a solid sphere has zero area and no incident flux.
        const bool is_center = is_boundary_face and std::fabs(cell_face.centroid.z) < 1.0e-12;

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
            else if (not is_center)
              psi = angle_set.PsiBoundary(cell_face.neighbor_id,
                                          direction_num,
                                          cell_local_id,
                                          f,
                                          fj,
                                          0,
                                          IsSurfaceSourceActive());
            if (not psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
          }
        }
      }

      const auto row_offset = static_cast<size_t>(direction_num) * num_moments_;
      const double* m2d_row = m2d_op.data() + row_offset;
      const double* d2m_row = d2m_op.data() + row_offset;

      for (size_t gsg = 0; gsg < gs_size && final_direction; ++gsg)
        for (size_t i = 0; i < cell_num_nodes; ++i)
          b[gsg](i) = psi_edge_[discretization_.MapDOFLocal(cell, i, unknown_manager_, 0, 0) + gsg];

      for (size_t gsg = 0; gsg < gs_size && not final_direction; ++gsg)
      {
        const double sigma_tg = sigma_t[gs_gi + gsg];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          double temp_src = 0.0;
          for (unsigned int m = 0; m < num_moments_; ++m)
            temp_src += m2d_row[m] * source_moments_[cell_transport_view.MapDOF(i, m, gs_gi + gsg)];
          source[i] = temp_src;
        }
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          double temp = 0.0;
          for (size_t j = 0; j < cell_num_nodes; ++j)
          {
            Atemp(i, j) = Amat(i, j) + M(i, j) * sigma_tg;
            temp += M(i, j) * source[j];
          }
          b[gsg](i) += temp;
        }
        GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes));
      }

      // Flux moments
      for (unsigned int m = 0; m < num_moments_; ++m)
      {
        const double wn_d2m = d2m_row[m];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const auto ir = cell_transport_view.MapDOF(i, m, gs_gi);
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            destination_phi_[ir + gsg] += wn_d2m * b[gsg](i);
        }
      }

      if (SaveAngularFluxEnabled())
      {
        double* cell_psi_data =
          &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const size_t imap =
            i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_psi_data[imap + gsg] = b[gsg](i);
        }
      }

      // Outgoing faces: downstream fluxes, reflected fluxes, and outflow tallies
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

          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            psi[gsg] = b[gsg](i);
        }
      }

      // Angular-edge flux for the next direction: psi_{n+1/2} = (psi_n - (1 - tau) psi_{n-1/2}) /
      // tau. The starting direction (tau = 1) initializes it.
      const double f0 = 1.0 / tau;
      const double f1 = f0 - 1.0;
      for (size_t i = 0; i < cell_num_nodes; ++i)
      {
        const auto ir = discretization_.MapDOFLocal(cell, i, unknown_manager_, 0, 0);
        for (size_t gsg = 0; gsg < gs_size; ++gsg)
          psi_edge_[ir + gsg] = f0 * b[gsg](i) - f1 * psi_edge_[ir + gsg];
      }
    }
  }
}

template <unsigned int NumNodes>
void
AAHSweepChunkSpherical::Sweep_FixedN(AngleSet& angle_set)
{
  constexpr std::size_t matrix_size = static_cast<std::size_t>(NumNodes) * NumNodes;
  constexpr auto idx = [](std::size_t i, std::size_t j) { return i * NumNodes + j; };
  const auto gs_size = groupset_.GetNumGroups();
  const auto gs_gi = groupset_.first_group;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();
  const auto& fac_diamond_difference = curvilinear_quadrature_.GetDiamondDifferenceFactor();
  const auto& fac_streaming_operator = curvilinear_quadrature_.GetStreamingOperatorFactor();

  std::vector<double> b(static_cast<std::size_t>(gs_size) * NumNodes);
  std::vector<double> face_mu_values;
  std::vector<double> sigma_block;
  sigma_block.reserve(group_block_size_);
  std::array<double, NumNodes> source{};
  std::vector<std::size_t> moment_dof_map(static_cast<std::size_t>(num_moments_) * NumNodes);

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const auto& as_angle_indices = angle_set.GetAngleIndices();
  for (size_t spls_index = 0; spls_index < spls.size(); ++spls_index)
  {
    const auto cell_local_id = spls[spls_index];
    const auto& cell = grid_->GetLocalCell(cell_local_id);
    const auto& cell_mapping = discretization_.GetCellMapping(cell);
    const auto& cell_transport_view = cell_transport_views_[cell_local_id];
    auto& cell_outflow_view = cell_outflow_views_[cell_local_id];
    const auto cell_num_faces = cell.faces.size();
    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    face_mu_values.resize(cell_num_faces);

    const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();
    const auto& G = unit_cell_matrices_[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = unit_cell_matrices_[cell_local_id].intS_shapeI_shapeJ;
    const auto& M_r = secondary_unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;

    std::array<double, matrix_size> mass_matrix{};
    std::array<double, matrix_size> angular_matrix{};
    std::array<std::uint64_t, NumNodes> edge_dof_map{};
    PRAGMA_UNROLL
    for (std::size_t i = 0; i < NumNodes; ++i)
    {
      edge_dof_map[i] = discretization_.MapDOFLocal(cell, i, unknown_manager_, 0, 0);
      PRAGMA_UNROLL
      for (std::size_t j = 0; j < NumNodes; ++j)
      {
        mass_matrix[idx(i, j)] = M(i, j);
        angular_matrix[idx(i, j)] = M_r(i, j);
      }
    }
    for (unsigned int m = 0; m < num_moments_; ++m)
    {
      PRAGMA_UNROLL
      for (std::size_t i = 0; i < NumNodes; ++i)
        moment_dof_map[static_cast<std::size_t>(m) * NumNodes + i] =
          cell_transport_view.MapDOF(i, m, gs_gi);
    }

    const auto ni_deploc_face_counter = deploc_face_counter;
    const auto ni_preloc_face_counter = preloc_face_counter;
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      const auto direction_num = as_angle_indices[as_ss_idx];
      const auto& omega = groupset_.quadrature->GetOmega(direction_num);
      const auto wt = groupset_.quadrature->GetWeight(direction_num);
      const double tau = fac_diamond_difference[direction_num];
      const double c_ang = fac_streaming_operator[direction_num];
      // The zero-weight final direction mu = +1 ends the recursion. Its flux is the angular-edge
      // flux psi_{N+1/2}; solving for it would be singular in a void cell at the center, whose
      // inflow face has zero area. It contributes nothing to the flux moments and is used only as
      // reflected inflow for the starting direction.
      const bool final_direction = wt == 0.0 and omega.z > 0.0;

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      // Angular-redistribution term: c_ang (psi_n - psi_{n-1/2}) / r
      std::fill(b.begin(), b.end(), 0.0);
      std::array<double, matrix_size> Amat{};
      PRAGMA_UNROLL
      for (std::size_t i = 0; i < NumNodes; ++i)
      {
        PRAGMA_UNROLL
        for (std::size_t j = 0; j < NumNodes; ++j)
        {
          const auto jr = edge_dof_map[j];
          const double streaming = c_ang * angular_matrix[idx(i, j)];
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            b[gsg * NumNodes + i] += streaming * psi_edge_[jr + gsg];
          Amat[idx(i, j)] = omega.Dot(G(i, j)) + streaming;
        }
      }

      for (size_t f = 0; f < cell_num_faces; ++f)
        face_mu_values[f] = omega.Dot(cell.faces[f].normal);

      // Incoming surface terms
      int in_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::INCOMING)
          continue;

        const auto& cell_face = cell.faces[f];
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary_face = not cell_face.has_neighbor;
        // The center of a solid sphere has zero area and no incident flux.
        const bool is_center = is_boundary_face and std::fabs(cell_face.centroid.z) < 1.0e-12;

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
            Amat[idx(i, j)] += mu_Nij;

            const double* psi = nullptr;
            if (is_local_face)
              psi = fluds.UpwindPsi(spls_index, in_face_counter, fj, 0, as_ss_idx);
            else if (not is_boundary_face)
              psi = fluds.NLUpwindPsi(preloc_face_counter, fj, 0, as_ss_idx);
            else if (not is_center)
              psi = angle_set.PsiBoundary(cell_face.neighbor_id,
                                          direction_num,
                                          cell_local_id,
                                          f,
                                          fj,
                                          0,
                                          IsSurfaceSourceActive());
            if (not psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg * NumNodes + i] += psi[gsg] * mu_Nij;
          }
        }
      }

      const auto row_offset = static_cast<size_t>(direction_num) * num_moments_;
      const double* m2d_row = m2d_op.data() + row_offset;
      const double* d2m_row = d2m_op.data() + row_offset;

      if (final_direction)
      {
        PRAGMA_UNROLL
        for (std::size_t i = 0; i < NumNodes; ++i)
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            b[gsg * NumNodes + i] = psi_edge_[edge_dof_map[i] + gsg];
      }

      for (unsigned int g0 = 0; g0 < gs_size and not final_direction; g0 += group_block_size_)
      {
        const auto g1 = std::min(g0 + group_block_size_, gs_size);
        const auto block_len = g1 - g0;
        sigma_block.resize(block_len);

        for (unsigned int gsg = g0; gsg < g1; ++gsg)
        {
          const auto rel = gsg - g0;
          sigma_block[rel] = sigma_t[gs_gi + gsg];

          source.fill(0.0);
          for (unsigned int m = 0; m < num_moments_; ++m)
          {
            const double w = m2d_row[m];
            PRAGMA_UNROLL
            for (std::size_t i = 0; i < NumNodes; ++i)
              source[i] +=
                w *
                source_moments_[moment_dof_map[static_cast<std::size_t>(m) * NumNodes + i] + gsg];
          }

          auto* __restrict bg = &b[static_cast<std::size_t>(gsg) * NumNodes];
          PRAGMA_UNROLL
          for (std::size_t i = 0; i < NumNodes; ++i)
          {
            double temp = 0.0;
            PRAGMA_UNROLL
            for (std::size_t j = 0; j < NumNodes; ++j)
              temp += mass_matrix[idx(i, j)] * source[j];
            bg[i] += temp;
          }
        }

        std::size_t k = 0;
        for (; k + Simd::size <= block_len; k += Simd::size)
          detail::SimdBatchSolve<NumNodes>(
            Amat.data(), mass_matrix.data(), &sigma_block[k], &b[(g0 + k) * NumNodes]);

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
      }

      // Flux moments
      for (unsigned int m = 0; m < num_moments_; ++m)
      {
        const double wn_d2m = d2m_row[m];
        PRAGMA_UNROLL
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
          const auto ir = moment_dof_map[static_cast<std::size_t>(m) * NumNodes + i];
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            destination_phi_[ir + gsg] += wn_d2m * b[gsg * NumNodes + i];
        }
      }

      if (SaveAngularFluxEnabled())
      {
        double* cell_psi_data =
          &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
        PRAGMA_UNROLL
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
          const size_t imap =
            i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_psi_data[imap + gsg] = b[gsg * NumNodes + i];
        }
      }

      // Outgoing faces: downstream fluxes, reflected fluxes, and outflow tallies
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
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_outflow_view.Add(
              f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg * NumNodes + i] * IntF_shapeI(i));

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
            psi[gsg] = b[gsg * NumNodes + i];
        }
      }

      // Angular-edge flux for the next direction: psi_{n+1/2} = (psi_n - (1 - tau) psi_{n-1/2}) /
      // tau. The starting direction (tau = 1) initializes it.
      const double f0 = 1.0 / tau;
      const double f1 = f0 - 1.0;
      PRAGMA_UNROLL
      for (std::size_t i = 0; i < NumNodes; ++i)
      {
        const auto ir = edge_dof_map[i];
        for (size_t gsg = 0; gsg < gs_size; ++gsg)
          psi_edge_[ir + gsg] = f0 * b[gsg * NumNodes + i] - f1 * psi_edge_[ir + gsg];
      }
    }
  }
}

template void AAHSweepChunkSpherical::Sweep_FixedN<2>(AngleSet&);

} // namespace opensn
