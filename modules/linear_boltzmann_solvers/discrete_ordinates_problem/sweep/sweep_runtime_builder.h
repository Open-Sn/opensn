// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace opensn
{

class GridFaceHistogram;
class LBSGroupset;
class MeshContinuum;
class SPDS;
class SpatialDiscretization;

using SweepOrderGroupingInfo = std::pair<UniqueSOGroupings, DirIDToSOMap>;

struct SweepRuntime
{
  std::map<std::shared_ptr<AngularQuadrature>, SweepOrderGroupingInfo>
    quadrature_unq_so_grouping_map;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::shared_ptr<SPDS>>>
    quadrature_spds_map;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::unique_ptr<FLUDSCommonData>>>
    quadrature_fluds_commondata_map;
};

SweepRuntime BuildSweepRuntime(const std::string& problem_name,
                               const std::vector<LBSGroupset>& groupsets,
                               const std::shared_ptr<MeshContinuum>& grid,
                               const std::string& sweep_type,
                               bool use_gpus,
                               const SpatialDiscretization& discretization,
                               const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

namespace detail
{

void BuildAAHGPUFludsCommonData(SweepRuntime& runtime,
                                const SpatialDiscretization& discretization,
                                const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

void BuildCBCGPUFludsCommonData(SweepRuntime& runtime,
                                const SpatialDiscretization& discretization,
                                const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

} // namespace detail

} // namespace opensn
