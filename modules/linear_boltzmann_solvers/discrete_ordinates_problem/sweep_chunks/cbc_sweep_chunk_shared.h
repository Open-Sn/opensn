// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_kernels.h"

namespace opensn
{

/// Reusable per-sweep-chunk context for binding angle-set and cell state.
struct CBCSweepChunkContext
{
  CBC_FLUDS* fluds = nullptr;
  CBC_AsynchronousCommunicator* async_comm = nullptr;
  size_t gs_size = 0;
  unsigned int gs_gi = 0;
  size_t num_angles_in_as = 0;
  unsigned int group_stride = 0;
  size_t group_angle_stride = 0;
  bool surface_source_active = false;
  const Cell* cell = nullptr;
  std::uint32_t cell_local_id = 0;
  const CellMapping* cell_mapping = nullptr;
  CellLBSView* cell_transport_view = nullptr;
  size_t cell_num_faces = 0;
  size_t cell_num_nodes = 0;
  const DenseMatrix<Vector3>* G = nullptr;
  const DenseMatrix<double>* M = nullptr;
  const std::vector<DenseMatrix<double>>* M_surf = nullptr;
  const std::vector<Vector<double>>* IntS_shapeI = nullptr;
  /// Reusable outgoing nonlocal face payload buffers.
  std::vector<CBCOutgoingFaceBuffer> outgoing_nonlocal_face_buffers;
  /// Outgoing nonlocal face payload lookup indexed by local face.
  std::vector<CBCOutgoingFaceBuffer*> outgoing_nonlocal_face_buffer_by_face;
  /// Number of outgoing nonlocal face payload buffers used by the current cell.
  size_t num_outgoing_nonlocal_face_buffers = 0;
  /// Reusable right-hand-side storage for fixed-node CBC kernels.
  std::vector<double> fixed_rhs_buffer;
  /// Reusable cross-section block storage for fixed-node CBC kernels.
  std::vector<double> fixed_sigma_block;
  /// Reusable moment/node DOF map for fixed-node CBC kernels.
  std::vector<size_t> fixed_moment_dof_map;
  /// Reusable angle-dependent face-normal dot products.
  std::vector<double> face_mu_values;
};

/// Bind angle-set state into a reusable CBC sweep-chunk context.
inline void
CBCBindAngleSetContext(CBCSweepChunkContext& ctx,
                       const LBSGroupset& groupset,
                       bool surface_source_active,
                       AngleSet& angle_set)
{
  ctx.fluds = &dynamic_cast<CBC_FLUDS&>(angle_set.GetFLUDS());
  ctx.async_comm = &dynamic_cast<CBC_AsynchronousCommunicator&>(*angle_set.GetCommunicator());
  ctx.gs_size = groupset.GetNumGroups();
  ctx.gs_gi = groupset.first_group;
  ctx.surface_source_active = surface_source_active;
  ctx.num_angles_in_as = angle_set.GetNumAngles();
  ctx.group_stride = angle_set.GetNumGroups();
  ctx.group_angle_stride = ctx.group_stride * ctx.num_angles_in_as;
}

/// Bind cell state into a reusable CBC sweep-chunk context.
inline void
CBCBindCellContext(CBCSweepChunkContext& ctx,
                   const SpatialDiscretization& discretization,
                   const std::vector<UnitCellMatrices>& unit_cell_matrices,
                   std::vector<CellLBSView>& cell_transport_views,
                   const Cell* cell_ptr)
{
  ctx.cell = cell_ptr;
  ctx.cell_local_id = cell_ptr->local_id;
  ctx.cell_mapping = &discretization.GetCellMapping(*ctx.cell);
  ctx.cell_transport_view = &cell_transport_views[ctx.cell->local_id];
  ctx.cell_num_faces = ctx.cell->faces.size();
  ctx.cell_num_nodes = ctx.cell_mapping->GetNumNodes();

  const auto& unit_mats = unit_cell_matrices[ctx.cell_local_id];
  ctx.G = &unit_mats.intV_shapeI_gradshapeJ;
  ctx.M = &unit_mats.intV_shapeI_shapeJ;
  ctx.M_surf = &unit_mats.intS_shapeI_shapeJ;
  ctx.IntS_shapeI = &unit_mats.intS_shapeI;
}

/// Create the bound data view consumed by CBC sweep kernels.
inline CBCSweepData
MakeCBCSweepData(const SpatialDiscretization& discretization,
                 const std::vector<double>& source_moments,
                 const LBSGroupset& groupset,
                 const BlockID2XSMap& xs,
                 unsigned int num_moments,
                 unsigned int max_num_cell_dofs,
                 bool save_angular_flux,
                 size_t groupset_angle_group_stride,
                 size_t groupset_group_stride,
                 std::vector<double>& destination_phi,
                 std::vector<double>& destination_psi,
                 bool include_rhs_time_term,
                 DiscreteOrdinatesProblem& problem,
                 const std::vector<double>* psi_old,
                 unsigned int group_block_size,
                 CBCSweepChunkContext& ctx)
{
  return CBCSweepData{discretization,
                      source_moments,
                      groupset,
                      xs,
                      num_moments,
                      max_num_cell_dofs,
                      save_angular_flux,
                      groupset_angle_group_stride,
                      groupset_group_stride,
                      destination_phi,
                      destination_psi,
                      ctx.surface_source_active,
                      include_rhs_time_term,
                      problem,
                      psi_old,
                      group_block_size,
                      *ctx.fluds,
                      *ctx.async_comm,
                      *ctx.cell,
                      ctx.cell_local_id,
                      *ctx.cell_mapping,
                      *ctx.cell_transport_view,
                      ctx.cell_num_faces,
                      ctx.cell_num_nodes,
                      ctx.gs_size,
                      ctx.gs_gi,
                      ctx.num_angles_in_as,
                      ctx.group_stride,
                      ctx.group_angle_stride,
                      *ctx.G,
                      *ctx.M,
                      *ctx.M_surf,
                      *ctx.IntS_shapeI,
                      ctx.outgoing_nonlocal_face_buffers,
                      ctx.outgoing_nonlocal_face_buffer_by_face,
                      ctx.num_outgoing_nonlocal_face_buffers,
                      ctx.fixed_rhs_buffer,
                      ctx.fixed_sigma_block,
                      ctx.fixed_moment_dof_map,
                      ctx.face_mu_values};
}

} // namespace opensn
