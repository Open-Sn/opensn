// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/discrete_ordinates_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/sweep_chunks/aah_sweep_chunk_rz.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <array>
#include <stdexcept>

namespace opensn
{

AAHSweepChunkRZ::AAHSweepChunkRZ(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
  : SweepChunk(problem.GetPhiNewLocal(),
               problem.GetPsiNewLocal()[groupset.id],
               problem.GetGrid(),
               problem.GetSpatialDiscretization(),
               problem.GetUnitCellMatrices(),
               problem.GetCellTransportViews(),
               problem.GetQMomentsLocal(),
               groupset,
               problem.GetBlockID2XSMap(),
               problem.GetNumMoments(),
               problem.GetMaxCellDOFCount(),
               problem.GetMinCellDOFCount()),
    secondary_unit_cell_matrices_(dynamic_cast<const DiscreteOrdinatesCurvilinearProblem&>(problem)
                                    .GetSecondaryUnitCellMatrices()),
    unknown_manager_(),
    psi_sweep_(),
    normal_vector_boundary_(),
    group_block_size_(ComputeGroupBlockSize(groupset.GetNumGroups())),
    sweep_impl_(&AAHSweepChunkRZ::Sweep_Generic)
{
  const auto curvilinear_product_quadrature =
    std::dynamic_pointer_cast<CurvilinearProductQuadrature>(groupset_.quadrature);

  if (curvilinear_product_quadrature == nullptr)
    throw std::invalid_argument("D_DO_RZ_SteadyState::SweepChunkPWL::SweepChunkPWL : "
                                "invalid angular quadrature");

  //  configure unknown manager for quantities that depend on polar level
  const size_t dir_map_size = curvilinear_product_quadrature->GetDirectionMap().size();
  for (size_t m = 0; m < dir_map_size; ++m)
    unknown_manager_.AddUnknown(UnknownType::VECTOR_N, groupset_.GetNumGroups());

  //  allocate storage for sweeping dependency
  const auto n_dof = discretization_.GetNumLocalDOFs(unknown_manager_);
  psi_sweep_.resize(n_dof);

  //  initialise mappings from direction linear index
  for (const auto& dir_set : curvilinear_product_quadrature->GetDirectionMap())
    for (const auto& dir_idx : dir_set.second)
      map_polar_level_.emplace(dir_idx, dir_set.first);

  //  set normal vector for symmetric boundary condition
  const auto d = (grid_->GetDimension() == 1) ? 2 : 0;
  normal_vector_boundary_ = Vector3(0.0, 0.0, 0.0);
  normal_vector_boundary_(d) = 1;

  if (min_num_cell_dofs_ == max_num_cell_dofs_)
  {
    switch (min_num_cell_dofs_)
    {
      case 2:
        sweep_impl_ = &AAHSweepChunkRZ::Sweep_FixedN<2>;
        break;
      case 3:
        sweep_impl_ = &AAHSweepChunkRZ::Sweep_FixedN<3>;
        break;
      case 4:
        sweep_impl_ = &AAHSweepChunkRZ::Sweep_FixedN<4>;
        break;
      case 5:
        sweep_impl_ = &AAHSweepChunkRZ::Sweep_FixedN<5>;
        break;
      case 6:
        sweep_impl_ = &AAHSweepChunkRZ::Sweep_FixedN<6>;
        break;
      case 7:
        sweep_impl_ = &AAHSweepChunkRZ::Sweep_FixedN<7>;
        break;
      case 8:
        sweep_impl_ = &AAHSweepChunkRZ::Sweep_FixedN<8>;
        break;
      default:
        break;
    }
  }
}

void
AAHSweepChunkRZ::Sweep(AngleSet& angle_set)
{
  (this->*sweep_impl_)(angle_set);
}

void
AAHSweepChunkRZ::Sweep_Generic(AngleSet& angle_set)
{
  auto gs_size = groupset_.GetNumGroups();
  auto gs_gi = groupset_.first_group;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();

  DenseMatrix<double> Amat(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Atemp(max_num_cell_dofs_, max_num_cell_dofs_);
  std::vector<Vector<double>> b(groupset_.GetNumGroups(), Vector<double>(max_num_cell_dofs_));
  std::vector<double> source(max_num_cell_dofs_);

  const auto curvilinear_product_quadrature =
    std::dynamic_pointer_cast<opensn::CurvilinearProductQuadrature>(groupset_.quadrature);

  // Loop over each cell
  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const size_t num_spls = spls.size();
  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    auto cell_local_id = spls[spls_index];
    const auto& cell = grid_->local_cells[cell_local_id];
    const auto& cell_mapping = discretization_.GetCellMapping(cell);
    auto& cell_transport_view = cell_transport_views_[cell_local_id];
    auto cell_num_faces = cell.faces.size();
    auto cell_num_nodes = cell_mapping.GetNumNodes();

    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    std::vector<double> face_mu_values(cell_num_faces);

    const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();

    // Get cell matrices
    const auto& G = unit_cell_matrices_[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = unit_cell_matrices_[cell_local_id].intS_shapeI_shapeJ;
    const auto& Maux = secondary_unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;

    // Loop over angles in set (as = angleset, ss = subset)
    const auto ni_deploc_face_counter = deploc_face_counter;
    const auto ni_preloc_face_counter = preloc_face_counter;
    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      auto direction_num = as_angle_indices[as_ss_idx];
      auto omega = groupset_.quadrature->omegas[direction_num];
      auto wt = groupset_.quadrature->weights[direction_num];

      const auto polar_level = map_polar_level_[direction_num];
      const auto fac_diamond_difference =
        curvilinear_product_quadrature->GetDiamondDifferenceFactor()[direction_num];
      const auto fac_streaming_operator =
        curvilinear_product_quadrature->GetStreamingOperatorFactor()[direction_num];

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      // Reset right-hand side
      for (size_t gsg = 0; gsg < gs_size; ++gsg)
        b[gsg] = Vector<double>(cell_num_nodes, 0.0);

      for (size_t i = 0; i < cell_num_nodes; ++i)
      {
        for (size_t j = 0; j < cell_num_nodes; ++j)
        {
          const auto jr =
            discretization_.MapDOFLocal(cell, j, unknown_manager_, polar_level, gs_gi);
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            b[gsg](i) += fac_streaming_operator * Maux(i, j) * psi_sweep_[jr + gsg];
        }
      }

      for (size_t i = 0; i < cell_num_nodes; ++i)
        for (size_t j = 0; j < cell_num_nodes; ++j)
          Amat(i, j) = omega.Dot(G(i, j)) + fac_streaming_operator * Maux(i, j);

      // Update face orientations
      for (size_t f = 0; f < cell_num_faces; ++f)
        face_mu_values[f] = omega.Dot(cell.faces[f].normal);

      // Surface integrals
      int in_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::INCOMING)
          continue;

        const auto& cell_face = cell.faces[f];
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary_face = not cell_face.has_neighbor;

        if (is_local_face)
          ++in_face_counter;
        else if (not is_boundary_face)
          ++preloc_face_counter;

        // IntSf_mu_psi_Mij_dA
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
            {
              //  Determine whether incoming direction is incident on the point
              //  of symmetry or on the axis of symmetry.
              //  N.B.: A face is considered to be on the point/axis of symmetry
              //  if all are true:
              //    1. The face normal is antiparallel to $\vec{e}_{d}$.
              //    2. All vertices of the face exhibit $v_{d} = 0$
              //       with $d = 2$ for 1D geometries and $d = 0$ for 2D geometries.
              //  Thanks to the verifications performed during initialisation,
              //  at this point it is necessary to confirm only the orientation.
              const bool incident_on_symmetric_boundary =
                (cell_face.normal.Dot(normal_vector_boundary_) < -0.999999);
              if (!incident_on_symmetric_boundary)
              {
                psi = angle_set.PsiBoundary(cell_face.neighbor_id,
                                            direction_num,
                                            cell_local_id,
                                            f,
                                            fj,
                                            gs_gi,
                                            IsSurfaceSourceActive());
              }
            }

            if (not psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
          } // for face node j
        } // for face node i
      } // for f

      const auto row_offset =
        static_cast<size_t>(direction_num) * static_cast<size_t>(num_moments_);
      const double* m2d_row = m2d_op.data() + row_offset;
      const double* d2m_row = d2m_op.data() + row_offset;

      // Looping over groups, assembling mass terms
      for (size_t gsg = 0; gsg < gs_size; ++gsg)
      {
        double sigma_tg = sigma_t[gs_gi + gsg];

        // Contribute source moments q = M_n^T * q_moms
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          double temp_src = 0.0;
          for (unsigned int m = 0; m < num_moments_; ++m)
          {
            const auto ir = cell_transport_view.MapDOF(i, m, gs_gi + gsg);
            temp_src += m2d_row[m] * source_moments_[ir];
          }
          source[i] = temp_src;
        }

        // Mass matrix and source
        // Atemp = Amat + sigma_tgr * M
        // b += M * q
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

        // Solve system
        GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes));
      } // for gsg

      // Update phi
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

      // Save angular flux during sweep
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

      // For outgoing, non-boundary faces, copy angular flux to fluds and
      // accumulate outflow
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
        } // for fi
      } // for face

      // Update sweeping dependency angular intensity for each polar level (incoming for next
      // interval)
      const auto f0 = 1 / fac_diamond_difference;
      const auto f1 = f0 - 1;
      for (size_t i = 0; i < cell_num_nodes; ++i)
      {
        const auto ir = discretization_.MapDOFLocal(cell, i, unknown_manager_, polar_level, gs_gi);
        for (size_t gsg = 0; gsg < gs_size; ++gsg)
          psi_sweep_[ir + gsg] = f0 * b[gsg](i) - f1 * psi_sweep_[ir + gsg];
      }
    } // for angleset/subset
  } // for cell
}

template <unsigned int NumNodes>
void
AAHSweepChunkRZ::Sweep_FixedN(AngleSet& angle_set)
{
  static_assert(NumNodes >= 2 and NumNodes <= 8);

  const auto gs_size = groupset_.GetNumGroups();
  const auto gs_gi = groupset_.first_group;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();

  constexpr size_t matrix_size = static_cast<size_t>(NumNodes) * NumNodes;
  constexpr auto idx = [](size_t i, size_t j) -> size_t { return i * NumNodes + j; };

  std::vector<double> b(static_cast<size_t>(gs_size) * NumNodes, 0.0);
  std::vector<double> sigma_block;
  sigma_block.reserve(group_block_size_);
  std::array<double, NumNodes> source{};
  std::vector<std::array<size_t, NumNodes>> moment_dof_map(num_moments_);

  const auto curvilinear_product_quadrature =
    std::dynamic_pointer_cast<CurvilinearProductQuadrature>(groupset_.quadrature);
  OpenSnLogicalErrorIf(curvilinear_product_quadrature == nullptr,
                       "AAHSweepChunkRZ requires a curvilinear product quadrature.");

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const size_t num_spls = spls.size();
  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    const auto cell_local_id = spls[spls_index];
    const auto& cell = grid_->local_cells[cell_local_id];
    const auto& cell_mapping = discretization_.GetCellMapping(cell);
    auto& cell_transport_view = cell_transport_views_[cell_local_id];
    const auto cell_num_faces = cell.faces.size();

    OpenSnInvalidArgumentIf(cell_mapping.GetNumNodes() != static_cast<size_t>(NumNodes),
                            "AAHSweepChunkRZ::Sweep_FixedN invoked for an incompatible cell "
                            "topology.");

    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    std::vector<double> face_mu_values(cell_num_faces);

    const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();

    const auto& G = unit_cell_matrices_[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = unit_cell_matrices_[cell_local_id].intS_shapeI_shapeJ;
    const auto& Maux = secondary_unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;

    std::array<double, matrix_size> mass_matrix{};
    std::array<double, matrix_size> aux_matrix{};
    PRAGMA_UNROLL
    for (size_t i = 0; i < NumNodes; ++i)
    {
      PRAGMA_UNROLL
      for (size_t j = 0; j < NumNodes; ++j)
      {
        mass_matrix[idx(i, j)] = M(i, j);
        aux_matrix[idx(i, j)] = Maux(i, j);
      }
    }

    for (unsigned int m = 0; m < num_moments_; ++m)
    {
      PRAGMA_UNROLL
      for (size_t i = 0; i < NumNodes; ++i)
        moment_dof_map[m][i] = cell_transport_view.MapDOF(i, m, gs_gi);
    }

    const auto ni_deploc_face_counter = deploc_face_counter;
    const auto ni_preloc_face_counter = preloc_face_counter;
    const auto& as_angle_indices = angle_set.GetAngleIndices();
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      const auto direction_num = as_angle_indices[as_ss_idx];
      const auto& omega = groupset_.quadrature->omegas[direction_num];
      const auto wt = groupset_.quadrature->weights[direction_num];

      const auto polar_level = map_polar_level_[direction_num];
      const auto fac_diamond_difference =
        curvilinear_product_quadrature->GetDiamondDifferenceFactor()[direction_num];
      const auto fac_streaming_operator =
        curvilinear_product_quadrature->GetStreamingOperatorFactor()[direction_num];

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      std::fill(b.begin(), b.end(), 0.0);

      std::array<double, matrix_size> Amat{};
      PRAGMA_UNROLL
      for (size_t i = 0; i < NumNodes; ++i)
      {
        PRAGMA_UNROLL
        for (size_t j = 0; j < NumNodes; ++j)
        {
          const auto jr =
            discretization_.MapDOFLocal(cell, j, unknown_manager_, polar_level, gs_gi);
          const double streaming = fac_streaming_operator * aux_matrix[idx(i, j)];
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            b[gsg * NumNodes + i] += streaming * psi_sweep_[jr + gsg];

          Amat[idx(i, j)] = omega.Dot(G(i, j)) + streaming;
        }
      }

      for (size_t f = 0; f < cell_num_faces; ++f)
        face_mu_values[f] = omega.Dot(cell.faces[f].normal);

      int in_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::INCOMING)
          continue;

        const auto& cell_face = cell.faces[f];
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
            Amat[idx(i, j)] += mu_Nij;

            const double* psi = nullptr;
            if (is_local_face)
              psi = fluds.UpwindPsi(spls_index, in_face_counter, fj, 0, as_ss_idx);
            else if (not is_boundary_face)
              psi = fluds.NLUpwindPsi(preloc_face_counter, fj, 0, as_ss_idx);
            else
            {
              const bool incident_on_symmetric_boundary =
                (cell_face.normal.Dot(normal_vector_boundary_) < -0.999999);
              if (not incident_on_symmetric_boundary)
              {
                psi = angle_set.PsiBoundary(cell_face.neighbor_id,
                                            direction_num,
                                            cell_local_id,
                                            f,
                                            fj,
                                            gs_gi,
                                            IsSurfaceSourceActive());
              }
            }

            if (psi != nullptr)
              for (size_t gsg = 0; gsg < gs_size; ++gsg)
                b[gsg * NumNodes + i] += psi[gsg] * mu_Nij;
          }
        }
      }

      const auto row_offset =
        static_cast<size_t>(direction_num) * static_cast<size_t>(num_moments_);
      const double* __restrict m2d_row = m2d_op.data() + row_offset;
      const double* __restrict d2m_row = d2m_op.data() + row_offset;

      for (unsigned int g0 = 0; g0 < gs_size; g0 += group_block_size_)
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
            for (size_t i = 0; i < NumNodes; ++i)
              source[i] += w * source_moments_[moment_dof_map[m][i] + gsg];
          }

          auto* __restrict bg = &b[static_cast<size_t>(gsg) * NumNodes];
          PRAGMA_UNROLL
          for (size_t i = 0; i < NumNodes; ++i)
          {
            double temp = 0.0;
            PRAGMA_UNROLL
            for (size_t j = 0; j < NumNodes; ++j)
              temp += mass_matrix[idx(i, j)] * source[j];
            bg[i] += temp;
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

        for (unsigned int gsg = g0; gsg < g1; ++gsg)
        {
          const auto* __restrict bg = &b[static_cast<size_t>(gsg) * NumNodes];
          for (unsigned int m = 0; m < num_moments_; ++m)
          {
            const double wn_d2m = d2m_row[m];
            PRAGMA_UNROLL
            for (size_t i = 0; i < NumNodes; ++i)
              destination_phi_[moment_dof_map[m][i] + gsg] += wn_d2m * bg[i];
          }
        }
      }

      if (SaveAngularFluxEnabled())
      {
        double* cell_psi_data =
          &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
        PRAGMA_UNROLL
        for (size_t i = 0; i < NumNodes; ++i)
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

        ++out_face_counter;
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
                f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg * NumNodes + i] * IntF_shapeI(i));
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
              psi[gsg] = b[gsg * NumNodes + i];
          }
        }
      }

      const auto f0 = 1.0 / fac_diamond_difference;
      const auto f1 = f0 - 1.0;
      PRAGMA_UNROLL
      for (size_t i = 0; i < NumNodes; ++i)
      {
        const auto ir = discretization_.MapDOFLocal(cell, i, unknown_manager_, polar_level, gs_gi);
        for (size_t gsg = 0; gsg < gs_size; ++gsg)
          psi_sweep_[ir + gsg] = f0 * b[gsg * NumNodes + i] - f1 * psi_sweep_[ir + gsg];
      }
    }
  }
}

template void AAHSweepChunkRZ::Sweep_FixedN<2>(AngleSet&);
template void AAHSweepChunkRZ::Sweep_FixedN<3>(AngleSet&);
template void AAHSweepChunkRZ::Sweep_FixedN<4>(AngleSet&);
template void AAHSweepChunkRZ::Sweep_FixedN<5>(AngleSet&);
template void AAHSweepChunkRZ::Sweep_FixedN<6>(AngleSet&);
template void AAHSweepChunkRZ::Sweep_FixedN<7>(AngleSet&);
template void AAHSweepChunkRZ::Sweep_FixedN<8>(AngleSet&);

} // namespace opensn
