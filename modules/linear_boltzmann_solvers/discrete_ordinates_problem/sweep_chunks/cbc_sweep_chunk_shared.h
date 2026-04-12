// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_kernels.h"

namespace opensn
{

/**
 * Cached per-cell and per-angleset context for CBC sweep chunks.
 *
 * Populated in two phases:
 * 1. BindAngleSet caches angleset-level data (FLUDS, group range, strides)
 * 2. BindCell caches cell-level data (geometry, transport views, unit cell matrices)
 */
struct CBCSweepChunkContext
{
  /// FLUDS for current angleset.
  CBC_FLUDS* fluds = nullptr;
  /// Number of groups in the current groupset.
  size_t gs_size = 0;
  /// First group index in the current groupset.
  unsigned int gs_gi = 0;
  /// Number of angles in the current angleset.
  size_t num_angles_in_as = 0;
  /// Per-angle group stride ( = num_groups).
  unsigned int group_stride = 0;
  /// Per-node angular stride ( = num_angles * num_groups).
  size_t group_angle_stride = 0;
  /// Whether the surface source BCs are active.
  bool surface_source_active = false;

  /// Current cell pointer
  const Cell* cell = nullptr;
  /// Local ID of the current cell.
  std::uint32_t cell_local_id = 0;
  /// Cell mapping for the current cell.
  const CellMapping* cell_mapping = nullptr;
  /// Transport view for the current cell.
  CellLBSView* cell_transport_view = nullptr;
  /// Number of faces on the current cell.
  size_t cell_num_faces = 0;
  /// Number of nodes in the current cell.
  size_t cell_num_nodes = 0;

  /// Volume integral matrix.
  const DenseMatrix<Vector3>* G = nullptr;
  /// Mass matrix.
  const DenseMatrix<double>* M = nullptr;
  /// Per-face surface mass matrices.
  const std::vector<DenseMatrix<double>>* M_surf = nullptr;
  /// Per-face surface integrals of shape functions.
  const std::vector<Vector<double>>* IntS_shapeI = nullptr;

  /// Cache angleset-level data (FLUDS, group range, strides).
  void BindAngleSet(const LBSGroupset& groupset, const bool has_surface_source, AngleSet& angle_set)
  {
    fluds = &dynamic_cast<CBC_FLUDS&>(angle_set.GetFLUDS());
    gs_size = groupset.GetNumGroups();
    gs_gi = groupset.first_group;
    surface_source_active = has_surface_source;
    num_angles_in_as = angle_set.GetNumAngles();
    group_stride = angle_set.GetNumGroups();
    group_angle_stride = group_stride * num_angles_in_as;
  }

  /// Cache cell-level data (geometry, transport views, unit cell matrices).
  void BindCell(const SpatialDiscretization& discretization,
                const std::vector<UnitCellMatrices>& unit_cell_matrices,
                std::vector<CellLBSView>& cell_transport_views,
                const Cell* cell_ptr)
  {
    cell = cell_ptr;
    cell_local_id = cell_ptr->local_id;
    cell_mapping = &discretization.GetCellMapping(*cell);
    cell_transport_view = &cell_transport_views[cell->local_id];
    cell_num_faces = cell->faces.size();
    cell_num_nodes = cell_mapping->GetNumNodes();

    const auto& unit_mats = unit_cell_matrices[cell_local_id];
    G = &unit_mats.intV_shapeI_gradshapeJ;
    M = &unit_mats.intV_shapeI_shapeJ;
    M_surf = &unit_mats.intS_shapeI_shapeJ;
    IntS_shapeI = &unit_mats.intS_shapeI;
  }
};

} // namespace opensn
