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
class Mesh;
class SPDS;
class SpatialDiscretization;

using SweepOrderGroupingInfo = std::pair<UniqueSOGroupings, DirIDToSOMap>;

struct SweepRuntime
{
  std::vector<std::shared_ptr<AngularQuadrature>> quadrature_order;
  std::map<std::shared_ptr<AngularQuadrature>, SweepOrderGroupingInfo>
    quadrature_unq_so_grouping_map;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::shared_ptr<SPDS>>>
    quadrature_spds_map;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::unique_ptr<FLUDSCommonData>>>
    quadrature_fluds_commondata_map;
};

SweepRuntime BuildSweepRuntime(const std::string& problem_name,
                               const std::vector<LBSGroupset>& groupsets,
                               const std::shared_ptr<Mesh>& grid,
                               const std::string& sweep_type,
                               bool use_gpus,
                               const SpatialDiscretization& discretization,
                               const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

namespace detail
{

/**
 * Group the directions of a quadrature by quadrant (`dimension` 2: signs of Omega_x, Omega_y) or
 * octant (`dimension` 3: signs of Omega_x, Omega_y, Omega_z). Bit d of a group's orthant index is
 * set when Omega_d is positive; empty groups are dropped, so the result has at most 2^dimension
 * groups ordered by orthant index. Directions keep their quadrature order within a group.
 *
 * Throws std::logic_error if `dimension` is not 2 or 3, or if a direction is not strictly inside
 * one quadrant/octant, i.e. |Omega_x| or |Omega_y| (or, in 3D, |Omega_z|) is at most
 * FACE_ORIENTATION_TOLERANCE.
 */
UniqueSOGroupings GroupDirectionsByOrthant(const std::string& problem_name,
                                           const AngularQuadrature& quadrature,
                                           unsigned int dimension);

/**
 * Direction from which the SPDS of a sweep ordering (a group of quadrature directions) is built.
 * The directions of a polar, quadrant or octant group give every face the same orientation on the
 * meshes those aggregation types allow, so their normalized mean direction gives the same sweep
 * ordering as any member, and its face weights stand for the whole group when cyclic dependencies
 * are broken. Other groups use their first direction.
 */
Vector3 SweepOrderingOmega(const AngularQuadrature& quadrature,
                           const DirIDs& so_grouping,
                           AngleAggregationType agg_type);

void BuildAAHGPUFludsCommonData(SweepRuntime& runtime,
                                const SpatialDiscretization& discretization,
                                const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

void BuildCBCGPUFludsCommonData(SweepRuntime& runtime,
                                const SpatialDiscretization& discretization,
                                const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

} // namespace detail

} // namespace opensn
