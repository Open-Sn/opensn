// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"

namespace opensn
{
namespace lbs
{
class LBSGroupset;

/** A sweep-chunk in point-symmetric and axial-symmetric
 *  curvilinear coordinates. */
class SweepChunkPwlrz : public SweepChunk
{
public:
  SweepChunkPwlrz(const MeshContinuum& grid,
                  const SpatialDiscretization& discretization_primary,
                  const std::vector<lbs::UnitCellMatrices>& unit_cell_matrices,
                  const std::vector<lbs::UnitCellMatrices>& secondary_unit_cell_matrices,
                  std::vector<lbs::CellLBSView>& cell_transport_views,
                  const std::vector<double>& densities,
                  std::vector<double>& destination_phi,
                  std::vector<double>& destination_psi,
                  const std::vector<double>& source_moments,
                  lbs::LBSGroupset& groupset,
                  const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                  int num_moments,
                  int max_num_cell_dofs);

  void Sweep(AngleSet& angle_set) override;

private:
  /** Secondary spatial discretization cell matrices */
  const std::vector<lbs::UnitCellMatrices>& secondary_unit_cell_matrices_;
  /** Unknown manager. */
  UnknownManager unknown_manager_;
  /** Sweeping dependency angular intensity (for each polar level). */
  std::vector<double> psi_sweep_;
  /** Mapping from direction linear index to direction polar level. */
  std::map<unsigned int, unsigned int> map_polar_level_;
  /** Normal vector to determine symmetric boundary condition. */
  Vector3 normal_vector_boundary_;
};

} // namespace lbs
} // namespace opensn
