// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"

namespace opensn
{

class AahSweepChunk : public SweepChunk
{
public:
  AahSweepChunk(const MeshContinuum& grid,
                const SpatialDiscretization& discretization,
                const std::vector<UnitCellMatrices>& unit_cell_matrices,
                std::vector<CellLBSView>& cell_transport_views,
                const std::vector<double>& densities,
                std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const std::vector<double>& source_moments,
                const LBSGroupset& groupset,
                const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                int num_moments,
                int max_num_cell_dofs);

  void Sweep(AngleSet& angle_set) override;
};

} // namespace opensn
