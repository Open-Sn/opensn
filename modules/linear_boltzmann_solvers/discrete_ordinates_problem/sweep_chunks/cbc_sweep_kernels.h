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

namespace opensn
{

/**
 * Aggregated sweep parameters for one cell in the CBC Generic/FixedN kernels.
 *
 * Bundles all data needed by the per-cell sweep kernel into a single struct
 * to avoid long parameter lists. Includes references to the spatial
 * discretization, groupset, FLUDS, cell geometry, and unit cell matrices.
 * Constructed once per cell by CBCSweepChunk::MakeSweepData.
 */
struct CBCSweepData
{
  /// Spatial discretization providing DOF mapping.
  const SpatialDiscretization& discretization;
  /// Source moment vector (indexed by DOF mapping).
  const std::vector<double>& source_moments;
  /// Groupset containing quadrature and group range.
  const LBSGroupset& groupset;
  /// Number of angular moments.
  unsigned int num_moments;
  /// Maximum number of DOFs (nodes) per cell in the mesh.
  unsigned int max_num_cell_dofs;
  /// Whether to store solved angular fluxes into destination_psi.
  bool save_angular_flux;
  /// Stride in the psi array: num_angles_in_quadrature * num_groups.
  size_t groupset_angle_group_stride;
  /// Stride in the psi array: num_groups.
  size_t groupset_group_stride;
  /// Output scalar flux moments vector (accumulated across angles).
  std::vector<double>& destination_phi;
  /// Output angular flux vector (written if save_angular_flux is true).
  std::vector<double>& destination_psi;
  /// Whether surface source boundary conditions are active.
  bool surface_source_active;
  /// Whether the RHS time-derivative term is included.
  bool include_rhs_time_term;
  /// Owning discrete ordinates problem (for time-step and theta access).
  DiscreteOrdinatesProblem& problem;
  /// Previous-time-step angular flux (nullptr for steady-state sweeps).
  const std::vector<double>* psi_old;
  /// Energy group block size for FixedN SIMD batch solve.
  unsigned int group_block_size;

  /// FLUDS providing local/nonlocal angular flux access.
  CBC_FLUDS& fluds;
  /// Current cell being swept.
  const Cell& cell;
  /// Local ID of the current cell.
  std::uint32_t cell_local_id;
  /// Cell mapping providing face-node maps and DOF counts.
  const CellMapping& cell_mapping;
  /// Transport view providing cross-section data and DOF mapping.
  CellLBSView& cell_transport_view;
  /// Number of faces on the current cell.
  size_t cell_num_faces;
  /// Number of nodes on the current cell.
  size_t cell_num_nodes;

  /// Number of energy groups in the groupset.
  size_t gs_size;
  /// First group index in the groupset.
  unsigned int gs_gi;
  /// Number of angles in the current angle set.
  size_t num_angles_in_as;
  /// Per-angle group stride (= num_groups).
  unsigned int group_stride;
  /// Per-node angular stride (= num_angles * num_groups).
  size_t group_angle_stride;

  /// Volume integral: \f$\int_V \nabla\phi_i \cdot \phi_j \, dV\f$.
  const DenseMatrix<Vector3>& G;
  /// Mass matrix: \f$\int_V \phi_i \phi_j \, dV\f$.
  const DenseMatrix<double>& M;
  /// Per-face surface mass matrices: \f$\int_S \phi_i \phi_j \, dS\f$.
  const std::vector<DenseMatrix<double>>& M_surf;
  /// Per-face surface integrals: \f$\int_S \phi_i \, dS\f$.
  const std::vector<Vector<double>>& IntS_shapeI;
};

/// Pre-resolved metadata for one incoming face of the current cell.
struct CBCIncomingFaceData
{
  /// Nodal mapping for local face access (nullptr for nonlocal/boundary).
  const FaceNodalMapping* face_nodal_mapping = nullptr;
  /// Base pointer for local/nonlocal incoming face psi; null only for boundary faces.
  double* psi_base = nullptr;
};

/// Pre-resolved metadata for one outgoing face of the current cell.
struct CBCOutgoingFaceData
{
  /// Whether the face is a reflecting boundary.
  bool is_reflecting_boundary_face = false;
  /// Base pointer for local outgoing face psi, when applicable.
  double* psi_base = nullptr;
  /// Nonlocal face info for MPI send staging; null for local/boundary faces.
  const CBC_FLUDSCommonData::OutgoingNonlocalFaceInfo* outgoing_nonlocal_face_info = nullptr;
};

/**
 * Reusable scratch buffers for the CBC Generic sweep kernel.
 *
 * Allocated once per sweep chunk and resized lazily via EnsureCapacity.
 * Avoids per-cell heap allocation in the hot path.
 */
struct CBCGenericSweepScratch
{
  /// Transport matrix: \f$A_{ij} = \hat\Omega \cdot G_{ij} + \text{face terms}\f$.
  DenseMatrix<double> Amat;
  /// Temporary copy of A with \f$\sigma_t M\f$ added, consumed by Gauss elimination.
  DenseMatrix<double> Atemp;
  /// Per-group RHS vectors.
  std::vector<Vector<double>> b;
  /// Per-node source assembly scratch.
  std::vector<double> source;
  /// Per-face dot product \f$\hat\Omega \cdot \hat n_f\f$.
  std::vector<double> face_mu_values;
  /// Per-group time-absorption coefficient \f$v_g^{-1} / (\theta \Delta t)\f$.
  std::vector<double> tau_gsg;
  /// Pre-resolved incoming face metadata (one per cell face).
  std::vector<CBCIncomingFaceData> incoming_face_data;
  /// Pre-resolved outgoing face metadata (one per cell face).
  std::vector<CBCOutgoingFaceData> outgoing_face_data;
  /// Pre-computed DOF indices: \c moment_dof_map[m * cell_num_nodes + i].
  std::vector<size_t> moment_dof_map;

  void
  EnsureCapacity(const size_t max_num_cell_dofs, const size_t gs_size, const size_t cell_num_faces)
  {
    if (Amat.Rows() != max_num_cell_dofs or Amat.Columns() != max_num_cell_dofs)
    {
      Amat = DenseMatrix<double>(max_num_cell_dofs, max_num_cell_dofs);
      Atemp = DenseMatrix<double>(max_num_cell_dofs, max_num_cell_dofs);
    }

    if (b.size() != gs_size)
      b.assign(gs_size, Vector<double>(max_num_cell_dofs));
    else
      for (auto& vec : b)
        if (vec.Rows() != max_num_cell_dofs)
          vec = Vector<double>(max_num_cell_dofs);

    if (source.size() != max_num_cell_dofs)
      source.assign(max_num_cell_dofs, 0.0);

    if (face_mu_values.size() != cell_num_faces)
      face_mu_values.assign(cell_num_faces, 0.0);

    if (incoming_face_data.size() != cell_num_faces)
      incoming_face_data.assign(cell_num_faces, CBCIncomingFaceData{});

    if (outgoing_face_data.size() != cell_num_faces)
      outgoing_face_data.assign(cell_num_faces, CBCOutgoingFaceData{});
  }
};

/**
 * Generic CBC sweep kernel for one cell, parameterized by time dependence.
 *
 * Assembles and solves the local transport system for all angles and groups
 * in the angle set, using dynamic-size matrices and Gauss elimination.
 * Used when the cell node count does not match a compile-time FixedN
 * specialization.
 *
 * \tparam time_dependent if true, include the time-derivative source term
 */
template <bool time_dependent>
inline void
CBC_Sweep_Generic(CBCSweepData& data, CBCGenericSweepScratch& scratch, AngleSet& angle_set)
{
  const auto& groupset = data.groupset;
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();
  scratch.EnsureCapacity(data.max_num_cell_dofs, data.gs_size, data.cell_num_faces);
  auto& Amat = scratch.Amat;
  auto& Atemp = scratch.Atemp;
  auto& b = scratch.b;
  auto& source = scratch.source;
  auto& face_mu_values = scratch.face_mu_values;

  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[data.cell_local_id];
  const auto& cell_xs = data.cell_transport_view.GetXS();
  const auto& sigma_t = cell_xs.GetSigmaTotal();

  scratch.tau_gsg.clear();
  if constexpr (time_dependent)
  {
    const auto& inv_velg = cell_xs.GetInverseVelocity();
    const double theta = data.problem.GetTheta();
    const double inv_theta = 1.0 / theta;
    const double dt = data.problem.GetTimeStep();
    const double inv_dt = 1.0 / dt;

    auto& tau_gsg = scratch.tau_gsg;
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
  auto& incoming_face_data = scratch.incoming_face_data;
  auto& outgoing_face_data = scratch.outgoing_face_data;
  for (size_t f = 0; f < data.cell_num_faces; ++f)
  {
    incoming_face_data[f] = CBCIncomingFaceData{};
    outgoing_face_data[f] = CBCOutgoingFaceData{};
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

  auto& moment_dof_map = scratch.moment_dof_map;
  moment_dof_map.resize(static_cast<size_t>(data.num_moments) * data.cell_num_nodes);
  for (unsigned int m = 0; m < data.num_moments; ++m)
    for (size_t i = 0; i < data.cell_num_nodes; ++i)
      moment_dof_map[static_cast<size_t>(m) * data.cell_num_nodes + i] =
        data.cell_transport_view.MapDOF(i, m, data.gs_gi);

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

    for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
      for (size_t i = 0; i < data.cell_num_nodes; ++i)
        b[gsg](i) = 0.0;

    for (size_t i = 0; i < data.cell_num_nodes; ++i)
      for (size_t j = 0; j < data.cell_num_nodes; ++j)
        Amat(i, j) = omega.Dot(data.G(i, j));

    for (size_t f = 0; f < data.cell_num_faces; ++f)
      face_mu_values[f] = omega.Dot(data.cell.faces[f].normal);

    for (size_t f = 0; f < data.cell_num_faces; ++f)
    {
      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const auto& face = data.cell.faces[f];
      const auto& face_data = incoming_face_data[f];
      const auto* face_nodal_mapping = face_data.face_nodal_mapping;

      const size_t num_face_nodes = data.cell_mapping.GetNumFaceNodes(f);
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = data.cell_mapping.MapFaceNode(f, fi);

        for (size_t fj = 0; fj < num_face_nodes; ++fj)
        {
          const int j = data.cell_mapping.MapFaceNode(f, fj);
          const double mu_Nij = -face_mu_values[f] * data.M_surf[f](i, j);
          Amat(i, j) += mu_Nij;

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

          if (psi != nullptr)
            for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
        }
      }
    }

    const auto dir_moment_offset =
      static_cast<std::size_t>(direction_num) * static_cast<std::size_t>(data.num_moments);
    const double* m2d_row = m2d_op.data() + dir_moment_offset;
    const double* d2m_row = d2m_op.data() + dir_moment_offset;

    for (unsigned int gsg = 0; gsg < data.gs_size; ++gsg)
    {
      double sigma_tg = sigma_t[data.gs_gi + gsg];
      if constexpr (time_dependent)
      {
        const auto& tau_gsg = scratch.tau_gsg;
        sigma_tg += tau_gsg[gsg];
      }

      for (size_t i = 0; i < data.cell_num_nodes; ++i)
      {
        double temp_src = 0.0;
        for (unsigned int m = 0; m < data.num_moments; ++m)
        {
          const auto ir = moment_dof_map[static_cast<size_t>(m) * data.cell_num_nodes + i] + gsg;
          temp_src += m2d_row[m] * data.source_moments[ir];
        }

        if constexpr (time_dependent)
        {
          const auto& tau_gsg = scratch.tau_gsg;
          const size_t imap =
            i * data.groupset_angle_group_stride + direction_num * data.groupset_group_stride;
          if (data.include_rhs_time_term and psi_old)
            temp_src += tau_gsg[gsg] * psi_old[imap + gsg];
        }

        source[i] = temp_src;
      }

      for (size_t i = 0; i < data.cell_num_nodes; ++i)
      {
        double temp = 0.0;
        for (size_t j = 0; j < data.cell_num_nodes; ++j)
        {
          const double Mij = data.M(i, j);
          Atemp(i, j) = Amat(i, j) + Mij * sigma_tg;
          temp += Mij * source[j];
        }
        b[gsg](i) += temp;
      }

      GaussElimination(Atemp, b[gsg], static_cast<int>(data.cell_num_nodes));
    }

    for (unsigned int m = 0; m < data.num_moments; ++m)
    {
      const auto wn_d2m = d2m_row[m];
      for (size_t i = 0; i < data.cell_num_nodes; ++i)
      {
        const auto ir = moment_dof_map[static_cast<size_t>(m) * data.cell_num_nodes + i];
        for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
          data.destination_phi[ir + gsg] += wn_d2m * b[gsg](i);
      }
    }

    if (data.save_angular_flux)
    {
      for (size_t i = 0; i < data.cell_num_nodes; ++i)
      {
        const size_t imap =
          i * data.groupset_angle_group_stride + direction_num * data.groupset_group_stride;

        for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
        {
          const double psi_sol = b[gsg](i);
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

      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = data.cell_mapping.MapFaceNode(f, fi);

        if (face_data.outgoing_nonlocal_face_info == nullptr and face_data.psi_base == nullptr)
        {
          for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
            data.cell_transport_view.AddOutflow(
              f, data.gs_gi + gsg, wt * face_mu_values[f] * b[gsg](i) * IntF_shapeI(i));
        }

        double* psi = nullptr;
        if (face_data.psi_base != nullptr)
          psi = face_data.psi_base + fi * data.group_angle_stride + as_ss_idx * data.group_stride;
        else if (face_data.outgoing_nonlocal_face_info != nullptr)
          psi = data.fluds.NLOutgoingPsi(psi_nonlocal_outgoing, fi, as_ss_idx);
        else if (is_reflecting_boundary_face)
          psi = angle_set.PsiReflected(face.neighbor_id, direction_num, data.cell_local_id, f, fi);

        if (psi != nullptr)
          for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
            psi[gsg] = b[gsg](i);
      }
    }
  }
}

/**
 * Fixed-node-count CBC sweep kernel with AVX/AVX512 SIMD batch solve.
 *
 * Specialized in cbc_avx_sweep_chunk.cc for compile-time-known node counts
 * (4, 8, etc.), enabling stack-allocated matrices, loop unrolling, and SIMD
 * batch Gauss elimination across multiple energy groups simultaneously.
 *
 * \tparam NumNodes compile-time number of cell nodes
 * \tparam time_dependent if true, include the time-derivative source term
 */
template <unsigned int NumNodes, bool time_dependent>
void CBC_Sweep_FixedN(CBCSweepData& data, AngleSet& angle_set);

} // namespace opensn
