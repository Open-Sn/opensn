// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/dense_matrix.h"
#include "framework/data_types/vector.h"
#include "framework/mesh/cell/cell.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include <algorithm>
#include <span>

namespace opensn
{

/// Staged outgoing nonlocal CBC face psi.
struct OutgoingFacePsi
{
  /// Downstream face slot.
  std::size_t face_slot = CBC_FLUDSCommonData::INVALID_FACE_SLOT;
  /// SPDS successor peer index.
  std::size_t peer_index = CBC_FLUDSCommonData::INVALID_PEER_INDEX;
  /// Downstream MPI rank for delayed face psi.
  int destination_location = -1;
  /// Whether this face psi is lagged by a broken inter-partition cycle.
  bool delayed = false;
  /// Face-node -> angleset subset -> group major psi values.
  std::vector<double> data;
};

/// Reusable host CBC sweep workspace.
struct CBCSweepWorkspace
{
  /// Staged nonlocal outgoing face psi.
  std::vector<OutgoingFacePsi> nonlocal_outgoing_psi;
  /// Face-indexed lookup for staged outgoing psi.
  std::vector<OutgoingFacePsi*> outgoing_psi_by_face;
  std::size_t num_nonlocal_outgoing_psi = 0;
  std::vector<double> rhs_values;
  std::vector<double> sigma_t_group_block;
  std::vector<std::size_t> moment_dof_map;
  std::vector<double> face_mu_values;
};

/// Prepare staged nonlocal outgoing face psi for one cell.
inline void
PrepareNonlocalOutgoingPsi(CBCSweepWorkspace& workspace,
                           CBC_FLUDS& fluds,
                           const Cell& cell,
                           std::uint32_t cell_local_id,
                           const CellMapping& cell_mapping,
                           const CellLBSView& cell_transport_view,
                           std::size_t group_angle_stride,
                           const std::vector<FaceOrientation>& face_orientations)
{
  auto& nonlocal_outgoing_psi = workspace.nonlocal_outgoing_psi;
  auto& outgoing_psi_by_face = workspace.outgoing_psi_by_face;
  const auto& common_data = fluds.GetCommonData();
  nonlocal_outgoing_psi.reserve(cell.faces.size());
  outgoing_psi_by_face.assign(cell.faces.size(), nullptr);
  workspace.num_nonlocal_outgoing_psi = 0;

  for (std::size_t f = 0; f < cell.faces.size(); ++f)
  {
    if (face_orientations[f] != FaceOrientation::OUTGOING)
      continue;

    const auto& face = cell.faces[f];
    if ((not face.has_neighbor) or cell_transport_view.IsFaceLocal(f))
      continue;

    const auto psi_index = workspace.num_nonlocal_outgoing_psi++;
    if (psi_index == nonlocal_outgoing_psi.size())
      nonlocal_outgoing_psi.emplace_back();

    auto& outgoing_psi = nonlocal_outgoing_psi[psi_index];
    outgoing_psi.face_slot =
      common_data.OutgoingFaceSlot(cell_local_id, static_cast<unsigned int>(f));
    outgoing_psi.delayed =
      common_data.IsDelayedNonlocalOutgoingFace(cell_local_id, static_cast<unsigned int>(f));
    outgoing_psi.destination_location =
      common_data.OutgoingFaceLocation(cell_local_id, static_cast<unsigned int>(f));
    outgoing_psi.peer_index =
      outgoing_psi.delayed
        ? CBC_FLUDSCommonData::INVALID_PEER_INDEX
        : common_data.OutgoingPeerIndex(cell_local_id, static_cast<unsigned int>(f));
    const auto psi_size = cell_mapping.GetNumFaceNodes(f) * group_angle_stride;
    if (outgoing_psi.data.size() != psi_size)
      outgoing_psi.data.resize(psi_size);
    outgoing_psi_by_face[f] = &outgoing_psi;
  }
}

/// Queue staged nonlocal outgoing face psi for asynchronous sends.
inline void
QueueNonlocalOutgoingPsi(CBCSweepWorkspace& workspace, CBC_AsynchronousCommunicator& async_comm)
{
  for (std::size_t i = 0; i < workspace.num_nonlocal_outgoing_psi; ++i)
  {
    const auto& outgoing_psi = workspace.nonlocal_outgoing_psi[i];
    async_comm.QueueDownwindMessage(
      outgoing_psi.delayed ? CBC_AsynchronousCommunicator::DownwindPsiType::DELAYED
                           : CBC_AsynchronousCommunicator::DownwindPsiType::NORMAL,
      outgoing_psi.delayed ? static_cast<std::size_t>(outgoing_psi.destination_location)
                           : outgoing_psi.peer_index,
      outgoing_psi.face_slot,
      std::span<const double>(outgoing_psi.data.data(), outgoing_psi.data.size()));
  }
}

/**
 * Sweep one host CBC cell using the generic dense-kernel path.
 * \tparam time_dependent Whether transient time terms are assembled.
 */
template <bool time_dependent, class SweepChunkT>
inline void
CBC_Sweep_Generic(SweepChunkT& sweep_chunk, AngleSet& angle_set)
{
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
  const std::size_t cell_num_nodes = cell_mapping.GetNumNodes();
  const auto& unit_mats = sweep_chunk.unit_cell_matrices_[cell_local_id];
  auto& fluds = *sweep_chunk.fluds_;
  const auto& common_data = fluds.GetCommonData();
  const auto& G = unit_mats.intV_shapeI_gradshapeJ;
  const auto& M = unit_mats.intV_shapeI_shapeJ;
  const auto& M_surf = unit_mats.intS_shapeI_shapeJ;
  const auto& IntS_shapeI = unit_mats.intS_shapeI;
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();

  DenseMatrix<double> Amat(sweep_chunk.max_num_cell_dofs_, sweep_chunk.max_num_cell_dofs_);
  DenseMatrix<double> Atemp(sweep_chunk.max_num_cell_dofs_, sweep_chunk.max_num_cell_dofs_);
  std::vector<Vector<double>> b(gs_size, Vector<double>(sweep_chunk.max_num_cell_dofs_));
  std::vector<double> source(sweep_chunk.max_num_cell_dofs_);
  std::vector<double> face_mu_values(cell_num_faces);

  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[cell_local_id];
  const auto& sigma_t = sweep_chunk.xs_.at(cell.block_id)->GetSigmaTotal();

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

    for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
      for (std::size_t i = 0; i < cell_num_nodes; ++i)
        b[gsg](i) = 0.0;

    for (std::size_t i = 0; i < cell_num_nodes; ++i)
      for (std::size_t j = 0; j < cell_num_nodes; ++j)
        Amat(i, j) = omega.Dot(G(i, j));

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

          if (psi != nullptr)
            for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
        }
      }
    }

    const auto dir_moment_offset =
      static_cast<std::size_t>(direction_num) * static_cast<std::size_t>(sweep_chunk.num_moments_);
    const double* m2d_row = m2d_op.data() + dir_moment_offset;
    const double* d2m_row = d2m_op.data() + dir_moment_offset;

    for (unsigned int gsg = 0; gsg < gs_size; ++gsg)
    {
      double sigma_tg = sigma_t[gs_gi + gsg];
      if constexpr (time_dependent)
        sigma_tg += tau_gsg[gsg];

      for (std::size_t i = 0; i < cell_num_nodes; ++i)
      {
        double temp_src = 0.0;
        for (unsigned int m = 0; m < sweep_chunk.num_moments_; ++m)
        {
          const auto ir = cell_transport_view.MapDOF(i, m, gs_gi + gsg);
          temp_src += m2d_row[m] * sweep_chunk.source_moments_[ir];
        }

        if constexpr (time_dependent)
        {
          const std::size_t imap = i * sweep_chunk.groupset_angle_group_stride_ +
                                   direction_num * sweep_chunk.groupset_group_stride_;
          if (sweep_chunk.include_rhs_time_term_ and psi_old)
            temp_src += tau_gsg[gsg] * psi_old[imap + gsg];
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

    for (unsigned int m = 0; m < sweep_chunk.num_moments_; ++m)
    {
      const auto wn_d2m = d2m_row[m];
      for (std::size_t i = 0; i < cell_num_nodes; ++i)
      {
        const auto ir = cell_transport_view.MapDOF(i, m, gs_gi);
        for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
          sweep_chunk.destination_phi_[ir + gsg] += wn_d2m * b[gsg](i);
      }
    }

    if (sweep_chunk.SaveAngularFluxEnabled())
    {
      double* psi_new = &sweep_chunk.destination_psi_[sweep_chunk.discretization_.MapDOFLocal(
        cell, 0, groupset.psi_uk_man_, 0, 0)];

      double theta = 1.0;
      double inv_theta = 1.0;
      if constexpr (time_dependent)
      {
        theta = sweep_chunk.problem_.GetTheta();
        inv_theta = 1.0 / theta;
      }

      for (std::size_t i = 0; i < cell_num_nodes; ++i)
      {
        const std::size_t imap = i * sweep_chunk.groupset_angle_group_stride_ +
                                 direction_num * sweep_chunk.groupset_group_stride_;

        for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
        {
          const double psi_sol = b[gsg](i);
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

      for (std::size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping.MapFaceNode(f, fi);

        if (is_boundary_face)
        {
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_outflow_view.Add(
              f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg](i) * IntF_shapeI(i));
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
          for (std::size_t gsg = 0; gsg < gs_size; ++gsg)
            psi[gsg] = b[gsg](i);
      }
    }
  }

  QueueNonlocalOutgoingPsi(sweep_chunk.workspace_, *sweep_chunk.async_comm_);
}

/**
 * Sweep one host CBC cell using a fixed-node-count dense-kernel path.
 * \tparam NumNodes Number of cell nodes.
 * \tparam time_dependent Whether transient time terms are assembled.
 */
template <unsigned int NumNodes, bool time_dependent, class SweepChunkT>
void CBC_Sweep_FixedN(SweepChunkT& sweep_chunk, AngleSet& angle_set);

} // namespace opensn
