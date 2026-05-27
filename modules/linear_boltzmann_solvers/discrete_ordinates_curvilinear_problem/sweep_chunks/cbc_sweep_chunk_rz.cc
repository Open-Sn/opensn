// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/sweep_chunks/cbc_sweep_chunk_rz.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/discrete_ordinates_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"
#include <algorithm>
#include <array>
#include <stdexcept>

namespace opensn
{
namespace
{

const CurvilinearProductQuadrature&
RequireCurvilinearProductQuadrature(const LBSGroupset& groupset)
{
  const auto* quadrature =
    dynamic_cast<const CurvilinearProductQuadrature*>(groupset.quadrature.get());
  if (quadrature == nullptr)
    throw std::invalid_argument("CBCSweepChunkRZ requires a curvilinear product quadrature.");
  return *quadrature;
}

} // namespace

CBCSweepChunkRZ::CBCSweepChunkRZ(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
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
    curvilinear_quadrature_(RequireCurvilinearProductQuadrature(groupset)),
    secondary_unit_cell_matrices_(dynamic_cast<const DiscreteOrdinatesCurvilinearProblem&>(problem)
                                    .GetSecondaryUnitCellMatrices()),
    unknown_manager_(),
    psi_sweep_(),
    direction_polar_level_(groupset.quadrature->GetNumAngles(), 0),
    normal_vector_boundary_(),
    group_block_size_(ComputeGroupBlockSize(groupset.GetNumGroups())),
    sweep_impl_(&CBCSweepChunkRZ::Sweep_Generic)
{
  const auto& direction_map = curvilinear_quadrature_.GetDirectionMap();
  const auto num_polar_levels = direction_map.size();
  for (std::size_t m = 0; m < num_polar_levels; ++m)
    unknown_manager_.AddUnknown(UnknownType::VECTOR_N, groupset_.GetNumGroups());

  psi_sweep_.assign(discretization_.GetNumLocalDOFs(unknown_manager_), 0.0);

  for (const auto& [polar_level, direction_indices] : direction_map)
  {
    for (const auto direction_index : direction_indices)
    {
      OpenSnLogicalErrorIf(direction_index >= direction_polar_level_.size(),
                           "CBCSweepChunkRZ received an invalid quadrature direction index.");
      direction_polar_level_[direction_index] = polar_level;
    }
  }

  const auto d = (grid_->GetDimension() == 1) ? 2 : 0;
  normal_vector_boundary_ = Vector3(0.0, 0.0, 0.0);
  normal_vector_boundary_(d) = 1.0;

  if (min_num_cell_dofs_ == max_num_cell_dofs_)
  {
    switch (min_num_cell_dofs_)
    {
      case 2:
        sweep_impl_ = &CBCSweepChunkRZ::Sweep_FixedN<2>;
        break;
      case 3:
        sweep_impl_ = &CBCSweepChunkRZ::Sweep_FixedN<3>;
        break;
      case 4:
        sweep_impl_ = &CBCSweepChunkRZ::Sweep_FixedN<4>;
        break;
      case 5:
        sweep_impl_ = &CBCSweepChunkRZ::Sweep_FixedN<5>;
        break;
      case 6:
        sweep_impl_ = &CBCSweepChunkRZ::Sweep_FixedN<6>;
        break;
      case 7:
        sweep_impl_ = &CBCSweepChunkRZ::Sweep_FixedN<7>;
        break;
      case 8:
        sweep_impl_ = &CBCSweepChunkRZ::Sweep_FixedN<8>;
        break;
      default:
        break;
    }
  }
}

void
CBCSweepChunkRZ::SetAngleSet(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunkRZ::SetAngleSet");

  fluds_ = &dynamic_cast<CBC_FLUDS&>(angle_set.GetFLUDS());
  async_comm_ = &dynamic_cast<CBC_AsynchronousCommunicator&>(*angle_set.GetCommunicator());
}

void
CBCSweepChunkRZ::Sweep(AngleSet& angle_set)
{
  (this->*sweep_impl_)(angle_set);
}

void
CBCSweepChunkRZ::Sweep_Generic(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunkRZ::Sweep");

  const auto& groupset = groupset_;
  const auto gs_size = groupset.GetNumGroups();
  const auto gs_gi = groupset.first_group;
  const auto num_angles_in_as = angle_set.GetNumAngles();
  const auto group_angle_stride = gs_size * num_angles_in_as;
  const auto& cell = *cell_;
  const auto cell_local_id = cell.local_id;
  const auto& cell_mapping = discretization_.GetCellMapping(cell);
  const auto& cell_transport_view = cell_transport_views_[cell_local_id];
  auto& cell_outflow_view = cell_outflow_views_[cell_local_id];
  const std::size_t cell_num_faces = cell.faces.size();
  const std::size_t cell_num_nodes = cell_mapping.GetNumNodes();
  const auto& unit_mats = unit_cell_matrices_[cell_local_id];
  auto& fluds = *fluds_;
  const auto& common_data = fluds.GetCommonData();
  const auto& G = unit_mats.intV_shapeI_gradshapeJ;
  const auto& M = unit_mats.intV_shapeI_shapeJ;
  const auto& M_surf = unit_mats.intS_shapeI_shapeJ;
  const auto& IntS_shapeI = unit_mats.intS_shapeI;
  const auto& Maux = secondary_unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();
  const auto& diamond_difference_factor = curvilinear_quadrature_.GetDiamondDifferenceFactor();
  const auto& streaming_operator_factor = curvilinear_quadrature_.GetStreamingOperatorFactor();
  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[cell_local_id];
  const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();

  DenseMatrix<double> Amat(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Atemp(max_num_cell_dofs_, max_num_cell_dofs_);
  std::vector<Vector<double>> b(gs_size, Vector<double>(max_num_cell_dofs_));
  std::vector<double> source(max_num_cell_dofs_);

  auto& face_mu_values = workspace_.face_mu_values;
  face_mu_values.assign(cell_num_faces, 0.0);
  PrepareNonlocalOutgoingPsi(workspace_,
                             fluds,
                             cell,
                             cell_local_id,
                             cell_mapping,
                             cell_transport_view,
                             group_angle_stride,
                             face_orientations);

  const auto& as_angle_indices = angle_set.GetAngleIndices();
  for (std::size_t as_ss_idx = 0; as_ss_idx < num_angles_in_as; ++as_ss_idx)
  {
    const auto direction_num = as_angle_indices[as_ss_idx];
    const auto& omega = groupset.quadrature->GetOmega(direction_num);
    const auto wt = groupset.quadrature->GetWeight(direction_num);

    const auto polar_level = direction_polar_level_[direction_num];
    const auto fac_diamond_difference = diamond_difference_factor[direction_num];
    const auto fac_streaming_operator = streaming_operator_factor[direction_num];

    for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
      for (std::size_t i = 0; i < cell_num_nodes; ++i)
        b[gsg](i) = 0.0;

    for (std::size_t i = 0; i < cell_num_nodes; ++i)
    {
      for (std::size_t j = 0; j < cell_num_nodes; ++j)
      {
        const auto jr = discretization_.MapDOFLocal(cell, j, unknown_manager_, polar_level, gs_gi);
        for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
          b[gsg](i) += fac_streaming_operator * Maux(i, j) * psi_sweep_[jr + gsg];
      }
    }

    for (std::size_t i = 0; i < cell_num_nodes; ++i)
      for (std::size_t j = 0; j < cell_num_nodes; ++j)
        Amat(i, j) = omega.Dot(G(i, j)) + fac_streaming_operator * Maux(i, j);

    for (std::size_t f = 0; f < cell_num_faces; ++f)
      face_mu_values[f] = omega.Dot(cell.faces[f].normal);

    for (std::size_t f = 0; f < cell_num_faces; ++f)
    {
      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const auto& face = cell.faces[f];
      const bool is_boundary_face = not face.has_neighbor;
      const bool is_local_face = not is_boundary_face and cell_transport_view.IsFaceLocal(f);
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

      const std::size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      for (std::size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping.MapFaceNode(f, fi);

        for (std::size_t fj = 0; fj < num_face_nodes; ++fj)
        {
          const int j = cell_mapping.MapFaceNode(f, fj);
          const double mu_Nij = -face_mu_values[f] * M_surf[f](i, j);
          Amat(i, j) += mu_Nij;

          const double* psi = nullptr;
          if (is_delayed_local_face)
          {
            psi = fluds.DelayedUpwindPsi(cell_local_id,
                                         static_cast<unsigned int>(f),
                                         face_nodal_mapping->face_node_mapping_[fj],
                                         as_ss_idx);
          }
          else if (is_delayed_nonlocal_face)
          {
            psi = fluds.DelayedNLUpwindPsi(
              delayed_nonlocal_face_info, face_nodal_mapping->face_node_mapping_[fj], as_ss_idx);
          }
          else if (is_local_face)
          {
            psi = fluds.UpwindPsi(*cell_transport_view.FaceNeighbor(f),
                                  face_nodal_mapping->cell_node_mapping_[fj],
                                  as_ss_idx);
          }
          else if (not is_boundary_face)
          {
            psi = fluds.NLUpwindPsi(
              incoming_nonlocal_slot, face_nodal_mapping->face_node_mapping_[fj], as_ss_idx);
          }
          else
          {
            const bool incident_on_symmetric_boundary =
              (face.normal.Dot(normal_vector_boundary_) < -0.999999);
            if (not incident_on_symmetric_boundary)
            {
              psi = angle_set.PsiBoundary(face.neighbor_id,
                                          direction_num,
                                          cell_local_id,
                                          static_cast<unsigned int>(f),
                                          static_cast<unsigned int>(fj),
                                          gs_gi,
                                          IsSurfaceSourceActive());
            }
          }

          if (psi != nullptr)
            for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
        }
      }
    }

    const auto dir_moment_offset =
      static_cast<std::size_t>(direction_num) * static_cast<std::size_t>(num_moments_);
    const double* m2d_row = m2d_op.data() + dir_moment_offset;
    const double* d2m_row = d2m_op.data() + dir_moment_offset;

    for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
    {
      const double sigma_tg = sigma_t[gs_gi + gsg];

      for (std::size_t i = 0; i < cell_num_nodes; ++i)
      {
        double temp_src = 0.0;
        for (unsigned int m = 0; m < num_moments_; ++m)
        {
          const auto ir = cell_transport_view.MapDOF(i, m, gs_gi + gsg);
          temp_src += m2d_row[m] * source_moments_[ir];
        }
        source[i] = temp_src;
      }

      for (std::size_t i = 0; i < cell_num_nodes; ++i)
      {
        double temp = 0.0;
        for (std::size_t j = 0; j < cell_num_nodes; ++j)
        {
          const double Mij = M(i, j);
          Atemp(i, j) = Amat(i, j) + Mij * sigma_tg;
          temp += Mij * source[j];
        }
        b[gsg](i) += temp;
      }

      GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes));
    }

    for (unsigned int m = 0; m < num_moments_; ++m)
    {
      const double wn_d2m = d2m_row[m];
      for (std::size_t i = 0; i < cell_num_nodes; ++i)
      {
        const auto ir = cell_transport_view.MapDOF(i, m, gs_gi);
        for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
          destination_phi_[ir + gsg] += wn_d2m * b[gsg](i);
      }
    }

    if (SaveAngularFluxEnabled())
    {
      double* cell_psi_data =
        &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)];

      for (std::size_t i = 0; i < cell_num_nodes; ++i)
      {
        const std::size_t imap =
          i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
        for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
          cell_psi_data[imap + gsg] = b[gsg](i);
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
      const auto& int_f_shape_i = IntS_shapeI[f];

      std::vector<double>* psi_nonlocal_outgoing = nullptr;
      if (not is_boundary_face and not is_local_face)
        psi_nonlocal_outgoing = &workspace_.outgoing_psi_by_face[f]->data;

      const std::size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      for (std::size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping.MapFaceNode(f, fi);

        if (is_boundary_face)
        {
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_outflow_view.Add(
              f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg](i) * int_f_shape_i(i));
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
        {
          psi = angle_set.PsiReflected(face.neighbor_id,
                                       direction_num,
                                       cell_local_id,
                                       static_cast<unsigned int>(f),
                                       static_cast<unsigned int>(fi));
        }

        if (psi != nullptr)
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            psi[gsg] = b[gsg](i);
      }
    }

    const auto f0 = 1.0 / fac_diamond_difference;
    const auto f1 = f0 - 1.0;
    for (std::size_t i = 0; i < cell_num_nodes; ++i)
    {
      const auto ir = discretization_.MapDOFLocal(cell, i, unknown_manager_, polar_level, gs_gi);
      for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
        psi_sweep_[ir + gsg] = f0 * b[gsg](i) - f1 * psi_sweep_[ir + gsg];
    }
  }

  QueueNonlocalOutgoingPsi(workspace_, *async_comm_);
}

template <unsigned int NumNodes>
void
CBCSweepChunkRZ::Sweep_FixedN(AngleSet& angle_set)
{
  static_assert(NumNodes >= 2 and NumNodes <= 8);

  const auto& groupset = groupset_;
  const auto gs_size = groupset.GetNumGroups();
  const auto gs_gi = groupset.first_group;
  const auto num_angles_in_as = angle_set.GetNumAngles();
  const auto group_angle_stride = gs_size * num_angles_in_as;
  const auto& cell = *cell_;
  const auto cell_local_id = cell.local_id;
  const auto& cell_mapping = discretization_.GetCellMapping(cell);
  const auto& cell_transport_view = cell_transport_views_[cell_local_id];
  auto& cell_outflow_view = cell_outflow_views_[cell_local_id];
  const std::size_t cell_num_faces = cell.faces.size();
  const auto& unit_mats = unit_cell_matrices_[cell_local_id];
  auto& fluds = *fluds_;
  const auto& common_data = fluds.GetCommonData();
  const auto& G = unit_mats.intV_shapeI_gradshapeJ;
  const auto& M = unit_mats.intV_shapeI_shapeJ;
  const auto& M_surf = unit_mats.intS_shapeI_shapeJ;
  const auto& IntS_shapeI = unit_mats.intS_shapeI;
  const auto& Maux = secondary_unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();
  const auto& diamond_difference_factor = curvilinear_quadrature_.GetDiamondDifferenceFactor();
  const auto& streaming_operator_factor = curvilinear_quadrature_.GetStreamingOperatorFactor();
  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[cell_local_id];
  const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();

  constexpr std::size_t matrix_size = static_cast<std::size_t>(NumNodes) * NumNodes;
  constexpr auto idx = [](std::size_t i, std::size_t j) -> std::size_t { return i * NumNodes + j; };

  std::array<double, matrix_size> mass_matrix{};
  std::array<double, matrix_size> aux_matrix{};
  PRAGMA_UNROLL
  for (std::size_t i = 0; i < NumNodes; ++i)
  {
    PRAGMA_UNROLL
    for (std::size_t j = 0; j < NumNodes; ++j)
    {
      mass_matrix[idx(i, j)] = M(i, j);
      aux_matrix[idx(i, j)] = Maux(i, j);
    }
  }

  auto& b = workspace_.rhs_values;
  b.resize(static_cast<std::size_t>(gs_size) * NumNodes);
  auto& sigma_block = workspace_.sigma_t_group_block;
  sigma_block.reserve(group_block_size_);
  std::array<double, NumNodes> source{};
  auto& moment_dof_map = workspace_.moment_dof_map;
  moment_dof_map.resize(static_cast<std::size_t>(num_moments_) * NumNodes);
  for (unsigned int m = 0; m < num_moments_; ++m)
  {
    PRAGMA_UNROLL
    for (std::size_t i = 0; i < NumNodes; ++i)
      moment_dof_map[static_cast<std::size_t>(m) * NumNodes + i] =
        cell_transport_view.MapDOF(i, m, gs_gi);
  }

  auto& face_mu_values = workspace_.face_mu_values;
  face_mu_values.assign(cell_num_faces, 0.0);
  PrepareNonlocalOutgoingPsi(workspace_,
                             fluds,
                             cell,
                             cell_local_id,
                             cell_mapping,
                             cell_transport_view,
                             group_angle_stride,
                             face_orientations);

  const auto& as_angle_indices = angle_set.GetAngleIndices();
  for (std::size_t as_ss_idx = 0; as_ss_idx < num_angles_in_as; ++as_ss_idx)
  {
    const auto direction_num = as_angle_indices[as_ss_idx];
    const auto& omega = groupset.quadrature->GetOmega(direction_num);
    const auto wt = groupset.quadrature->GetWeight(direction_num);

    const auto polar_level = direction_polar_level_[direction_num];
    const auto fac_diamond_difference = diamond_difference_factor[direction_num];
    const auto fac_streaming_operator = streaming_operator_factor[direction_num];

    std::fill(b.begin(), b.end(), 0.0);

    std::array<double, matrix_size> Amat{};
    PRAGMA_UNROLL
    for (std::size_t i = 0; i < NumNodes; ++i)
    {
      PRAGMA_UNROLL
      for (std::size_t j = 0; j < NumNodes; ++j)
      {
        const auto jr = discretization_.MapDOFLocal(cell, j, unknown_manager_, polar_level, gs_gi);
        const double streaming = fac_streaming_operator * aux_matrix[idx(i, j)];
        for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
          b[gsg * NumNodes + i] += streaming * psi_sweep_[jr + gsg];

        Amat[idx(i, j)] = omega.Dot(G(i, j)) + streaming;
      }
    }

    for (std::size_t f = 0; f < cell_num_faces; ++f)
      face_mu_values[f] = omega.Dot(cell.faces[f].normal);

    for (std::size_t f = 0; f < cell_num_faces; ++f)
    {
      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const auto& face = cell.faces[f];
      const bool is_boundary_face = not face.has_neighbor;
      const bool is_local_face = not is_boundary_face and cell_transport_view.IsFaceLocal(f);
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

      const std::size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      for (std::size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping.MapFaceNode(f, fi);
        for (std::size_t fj = 0; fj < num_face_nodes; ++fj)
        {
          const int j = cell_mapping.MapFaceNode(f, fj);
          const double mu_Nij = -face_mu_values[f] * M_surf[f](i, j);
          Amat[idx(i, j)] += mu_Nij;

          const double* psi = nullptr;
          if (is_delayed_local_face)
          {
            psi = fluds.DelayedUpwindPsi(cell_local_id,
                                         static_cast<unsigned int>(f),
                                         face_nodal_mapping->face_node_mapping_[fj],
                                         as_ss_idx);
          }
          else if (is_delayed_nonlocal_face)
          {
            psi = fluds.DelayedNLUpwindPsi(
              delayed_nonlocal_face_info, face_nodal_mapping->face_node_mapping_[fj], as_ss_idx);
          }
          else if (is_local_face)
          {
            psi = fluds.UpwindPsi(*cell_transport_view.FaceNeighbor(f),
                                  face_nodal_mapping->cell_node_mapping_[fj],
                                  as_ss_idx);
          }
          else if (not is_boundary_face)
          {
            psi = fluds.NLUpwindPsi(
              incoming_nonlocal_slot, face_nodal_mapping->face_node_mapping_[fj], as_ss_idx);
          }
          else
          {
            const bool incident_on_symmetric_boundary =
              (face.normal.Dot(normal_vector_boundary_) < -0.999999);
            if (not incident_on_symmetric_boundary)
            {
              psi = angle_set.PsiBoundary(face.neighbor_id,
                                          direction_num,
                                          cell_local_id,
                                          static_cast<unsigned int>(f),
                                          static_cast<unsigned int>(fj),
                                          gs_gi,
                                          IsSurfaceSourceActive());
            }
          }

          if (psi != nullptr)
            for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg * NumNodes + i] += psi[gsg] * mu_Nij;
        }
      }
    }

    const auto dir_moment_offset =
      static_cast<std::size_t>(direction_num) * static_cast<std::size_t>(num_moments_);
    const double* __restrict m2d_row = m2d_op.data() + dir_moment_offset;
    const double* __restrict d2m_row = d2m_op.data() + dir_moment_offset;

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
          const auto moment_offset = static_cast<std::size_t>(m) * NumNodes;
          PRAGMA_UNROLL
          for (std::size_t i = 0; i < NumNodes; ++i)
            source[i] += w * source_moments_[moment_dof_map[moment_offset + i] + gsg];
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

      for (unsigned int gsg = g0; gsg < g1; ++gsg)
      {
        const auto* __restrict bg = &b[static_cast<std::size_t>(gsg) * NumNodes];
        for (unsigned int m = 0; m < num_moments_; ++m)
        {
          const double wn_d2m = d2m_row[m];
          const auto moment_offset = static_cast<std::size_t>(m) * NumNodes;
          PRAGMA_UNROLL
          for (std::size_t i = 0; i < NumNodes; ++i)
            destination_phi_[moment_dof_map[moment_offset + i] + gsg] += wn_d2m * bg[i];
        }
      }
    }

    if (SaveAngularFluxEnabled())
    {
      double* cell_psi_data =
        &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)];

      PRAGMA_UNROLL
      for (std::size_t i = 0; i < NumNodes; ++i)
      {
        const std::size_t imap =
          i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
        for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
          cell_psi_data[imap + gsg] = b[gsg * NumNodes + i];
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
      const auto& int_f_shape_i = IntS_shapeI[f];

      std::vector<double>* psi_nonlocal_outgoing = nullptr;
      if (not is_boundary_face and not is_local_face)
        psi_nonlocal_outgoing = &workspace_.outgoing_psi_by_face[f]->data;

      const std::size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      for (std::size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping.MapFaceNode(f, fi);

        if (is_boundary_face)
        {
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_outflow_view.Add(
              f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg * NumNodes + i] * int_f_shape_i(i));
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
        {
          psi = angle_set.PsiReflected(face.neighbor_id,
                                       direction_num,
                                       cell_local_id,
                                       static_cast<unsigned int>(f),
                                       static_cast<unsigned int>(fi));
        }

        if (psi != nullptr)
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            psi[gsg] = b[gsg * NumNodes + i];
      }
    }

    const auto f0 = 1.0 / fac_diamond_difference;
    const auto f1 = f0 - 1.0;
    PRAGMA_UNROLL
    for (std::size_t i = 0; i < NumNodes; ++i)
    {
      const auto ir = discretization_.MapDOFLocal(cell, i, unknown_manager_, polar_level, gs_gi);
      for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
        psi_sweep_[ir + gsg] = f0 * b[gsg * NumNodes + i] - f1 * psi_sweep_[ir + gsg];
    }
  }

  QueueNonlocalOutgoingPsi(workspace_, *async_comm_);
}

template void CBCSweepChunkRZ::Sweep_FixedN<2>(AngleSet&);
template void CBCSweepChunkRZ::Sweep_FixedN<3>(AngleSet&);
template void CBCSweepChunkRZ::Sweep_FixedN<4>(AngleSet&);
template void CBCSweepChunkRZ::Sweep_FixedN<5>(AngleSet&);
template void CBCSweepChunkRZ::Sweep_FixedN<6>(AngleSet&);
template void CBCSweepChunkRZ::Sweep_FixedN<7>(AngleSet&);
template void CBCSweepChunkRZ::Sweep_FixedN<8>(AngleSet&);

} // namespace opensn
