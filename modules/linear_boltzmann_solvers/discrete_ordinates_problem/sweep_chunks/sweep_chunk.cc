// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/utils/error.h"
#include <algorithm>

namespace opensn
{
void
SweepChunk::Sweep(AngleSet& angle_set)
{
  OpenSnLogicalError("SweepChunk::Sweep must be overridden by derived sweep chunk.");
}

void
SweepChunk::SetAngleSet(AngleSet& angle_set)
{
  OpenSnLogicalError("SweepChunk::SetAngleSet is not implemented for this sweep chunk.");
}

void
SweepChunk::SetCell(Cell const* cell_ptr)
{
  cell_ = cell_ptr;
}

void
SweepChunk::ZeroDestinationPhi()
{
  const auto gsi = groupset_.first_group;
  const auto gss = groupset_.GetNumGroups();

  for (const auto& cell : grid_->GetLocalCells())
  {
    const auto& transport_view = cell_transport_views_[cell->local_id];
    const auto num_nodes = static_cast<size_t>(transport_view.GetNumNodes());
    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (unsigned int m = 0; m < num_moments_; ++m)
      {
        const auto mapping = transport_view.MapDOF(i, m, gsi);
        for (unsigned int g = 0; g < gss; ++g)
        {
          destination_phi_[mapping + g] = 0.0;
        } // for g
      } // for moment
    } // for dof

  } // for cell
}

void
SweepChunk::SetDestinationsFromAngularFlux(const AngleAggregation& angle_agg)
{
  ComputeDestinationsFromAngularFlux(angle_agg, destination_psi_);
}

void
SweepChunk::ComputeDestinationsFromAngularFlux(const AngleAggregation& angle_agg,
                                               const std::vector<double>& psi)
{
  OpenSnLogicalErrorIf(psi.empty(),
                       "Rebuilding sweep destinations requires stored angular fluxes.");

  const auto& quadrature = *groupset_.quadrature;
  const auto& d2m_op = quadrature.GetDiscreteToMomentOperator();
  const auto gs_gi = groupset_.first_group;
  const auto gs_size = groupset_.GetNumGroups();

  for (const auto& cell : grid_->GetLocalCells())
  {
    const auto cell_local_id = cell->local_id;
    const auto& transport_view = cell_transport_views_[cell_local_id];
    const auto& cell_mapping = discretization_.GetCellMapping(*cell);
    const auto& int_s_shape_i = unit_cell_matrices_[cell_local_id].intS_shapeI;
    auto& outflow_view = cell_outflow_views_[cell_local_id];
    const auto num_nodes = static_cast<size_t>(transport_view.GetNumNodes());
    const auto num_faces = cell->faces.size();
    const double* cell_psi =
      &psi[discretization_.MapDOFLocal(*cell, 0, groupset_.psi_uk_man_, 0, 0)];

    for (size_t i = 0; i < num_nodes; ++i)
      for (unsigned int m = 0; m < num_moments_; ++m)
        std::fill_n(&destination_phi_[transport_view.MapDOF(i, m, gs_gi)], gs_size, 0.0);
    for (size_t f = 0; f < num_faces; ++f)
      for (unsigned int gsg = 0; gsg < gs_size; ++gsg)
        outflow_view.Zero(f, gs_gi + gsg);

    // Same terms as the sweep chunks, summed in a fixed order.
    for (const auto& angle_set : angle_agg)
    {
      const auto& face_orientations = angle_set->GetSPDS().GetCellFaceOrientations()[cell_local_id];
      for (const auto direction_num : angle_set->GetAngleIndices())
      {
        const double* d2m_row = d2m_op.data() + direction_num * static_cast<size_t>(num_moments_);
        for (size_t i = 0; i < num_nodes; ++i)
        {
          const double* psi_i =
            cell_psi + i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
          for (unsigned int m = 0; m < num_moments_; ++m)
          {
            const auto wn_d2m = d2m_row[m];
            double* phi = &destination_phi_[transport_view.MapDOF(i, m, gs_gi)];
            for (unsigned int gsg = 0; gsg < gs_size; ++gsg)
              phi[gsg] += wn_d2m * psi_i[gsg];
          }
        }

        const auto& omega = quadrature.GetOmega(direction_num);
        const auto wt = quadrature.GetWeight(direction_num);
        for (size_t f = 0; f < num_faces; ++f)
        {
          if (face_orientations[f] != FaceOrientation::OUTGOING)
            continue;
          const double mu = omega.Dot(cell->faces[f].normal);
          const auto& int_f_shape_i = int_s_shape_i[f];
          for (size_t fi = 0; fi < cell_mapping.GetNumFaceNodes(f); ++fi)
          {
            const int i = cell_mapping.MapFaceNode(f, fi);
            const double* psi_i =
              cell_psi + i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
            for (unsigned int gsg = 0; gsg < gs_size; ++gsg)
              outflow_view.Add(f, gs_gi + gsg, wt * mu * psi_i[gsg] * int_f_shape_i(i));
          }
        }
      }
    }
  } // for cell
}

void
SweepChunk::ComputeTimeDependentDestinationsFromAngularFlux(const AngleAggregation& angle_agg,
                                                            const std::vector<double>& psi_old,
                                                            const double theta)
{
  OpenSnLogicalErrorIf(psi_old.size() != destination_psi_.size(),
                       "Time-dependent sweep destinations require the previous angular flux.");
  psi_theta_.resize(destination_psi_.size());
  for (size_t k = 0; k < psi_theta_.size(); ++k)
    psi_theta_[k] = theta * destination_psi_[k] + (1.0 - theta) * psi_old[k];
  ComputeDestinationsFromAngularFlux(angle_agg, psi_theta_);
}

} // namespace opensn
