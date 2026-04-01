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

struct CBCSweepData
{
  const SpatialDiscretization& discretization;
  const std::vector<double>& source_moments;
  const LBSGroupset& groupset;
  const BlockID2XSMap& xs;
  unsigned int num_moments;
  unsigned int max_num_cell_dofs;
  bool save_angular_flux;
  size_t groupset_angle_group_stride;
  size_t groupset_group_stride;
  std::vector<double>& destination_phi;
  std::vector<double>& destination_psi;
  bool surface_source_active;
  bool include_rhs_time_term;
  DiscreteOrdinatesProblem& problem;
  const std::vector<double>* psi_old;
  unsigned int group_block_size;

  CBC_FLUDS& fluds;
  const Cell& cell;
  std::uint32_t cell_local_id;
  const CellMapping& cell_mapping;
  CellLBSView& cell_transport_view;
  size_t cell_num_faces;
  size_t cell_num_nodes;

  size_t gs_size;
  unsigned int gs_gi;
  size_t num_angles_in_as;
  unsigned int group_stride;
  size_t group_angle_stride;

  const DenseMatrix<Vector3>& G;
  const DenseMatrix<double>& M;
  const std::vector<DenseMatrix<double>>& M_surf;
  const std::vector<Vector<double>>& IntS_shapeI;
};

template <bool time_dependent>
inline void
CBC_Sweep_Generic(CBCSweepData& data, AngleSet& angle_set)
{
  const auto& groupset = data.groupset;
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();

  DenseMatrix<double> Amat(data.max_num_cell_dofs, data.max_num_cell_dofs);
  DenseMatrix<double> Atemp(data.max_num_cell_dofs, data.max_num_cell_dofs);
  std::vector<Vector<double>> b(data.gs_size, Vector<double>(data.max_num_cell_dofs));
  std::vector<double> source(data.max_num_cell_dofs);
  std::vector<double> face_mu_values(data.cell_num_faces);

  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[data.cell_local_id];
  const auto& sigma_t = data.xs.at(data.cell.block_id)->GetSigmaTotal();

  std::vector<double> tau_gsg;
  if constexpr (time_dependent)
  {
    const auto& inv_velg = data.xs.at(data.cell.block_id)->GetInverseVelocity();
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
      const bool is_local_face = data.cell_transport_view.IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      const auto* face_nodal_mapping =
        &data.fluds.GetCommonData().GetFaceNodalMapping(data.cell_local_id, f);

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

          if (is_local_face)
            psi = data.fluds.UpwindPsi(*data.cell_transport_view.FaceNeighbor(f),
                                       face_nodal_mapping->cell_node_mapping_[fj],
                                       as_ss_idx);
          else if (not is_boundary_face)
            psi = data.fluds.NLUpwindPsi(
              data.cell.global_id, f, face_nodal_mapping->face_node_mapping_[fj], as_ss_idx);
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
        sigma_tg += tau_gsg[gsg];

      for (size_t i = 0; i < data.cell_num_nodes; ++i)
      {
        double temp_src = 0.0;
        for (unsigned int m = 0; m < data.num_moments; ++m)
        {
          const auto ir = data.cell_transport_view.MapDOF(i, m, data.gs_gi + gsg);
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
        const auto ir = data.cell_transport_view.MapDOF(i, m, data.gs_gi);
        for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
          data.destination_phi[ir + gsg] += wn_d2m * b[gsg](i);
      }
    }

    if (data.save_angular_flux)
    {
      double* psi_new = &data.destination_psi[data.discretization.MapDOFLocal(
        data.cell, 0, groupset.psi_uk_man_, 0, 0)];

      double theta = 1.0;
      double inv_theta = 1.0;
      if constexpr (time_dependent)
      {
        theta = data.problem.GetTheta();
        inv_theta = 1.0 / theta;
      }

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
            psi_new[imap + gsg] = inv_theta * (psi_sol + (theta - 1.0) * psi_old_val);
          }
          else
            psi_new[imap + gsg] = psi_sol;
        }
      }
    }

    for (size_t f = 0; f < data.cell_num_faces; ++f)
    {
      if (face_orientations[f] != FaceOrientation::OUTGOING)
        continue;

      const auto& face = data.cell.faces[f];
      const bool is_local_face = data.cell_transport_view.IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      const bool is_reflecting_boundary_face =
        (is_boundary_face and angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting());
      const auto& IntF_shapeI = data.IntS_shapeI[f];

      const int locality = data.cell_transport_view.FaceLocality(f);
      const size_t num_face_nodes = data.cell_mapping.GetNumFaceNodes(f);
      const auto& face_nodal_mapping =
        data.fluds.GetCommonData().GetFaceNodalMapping(data.cell_local_id, f);
      std::vector<double>* psi_nonlocal_outgoing = nullptr;

      if (not is_boundary_face and not is_local_face)
      {
        auto* async_comm = dynamic_cast<CBC_AsynchronousCommunicator*>(angle_set.GetCommunicator());
        const size_t data_size_for_msg = num_face_nodes * data.group_angle_stride;
        psi_nonlocal_outgoing =
          &async_comm->InitGetDownwindMessageData(locality,
                                                  face.neighbor_id,
                                                  face_nodal_mapping.associated_face_,
                                                  angle_set.GetID(),
                                                  data_size_for_msg);
      }

      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = data.cell_mapping.MapFaceNode(f, fi);

        if (is_boundary_face)
        {
          for (size_t gsg = 0; gsg < data.gs_size; ++gsg)
            data.cell_transport_view.AddOutflow(
              f, data.gs_gi + gsg, wt * face_mu_values[f] * b[gsg](i) * IntF_shapeI(i));
        }

        double* psi = nullptr;
        if (is_local_face)
          psi = data.fluds.OutgoingPsi(data.cell, i, as_ss_idx);
        else if (not is_boundary_face)
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

template <unsigned int NumNodes, bool time_dependent>
void CBC_Sweep_FixedN(CBCSweepData& data, AngleSet& angle_set);

} // namespace opensn
