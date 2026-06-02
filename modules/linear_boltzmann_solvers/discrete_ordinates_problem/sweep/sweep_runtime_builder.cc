// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep_runtime_builder.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/logging/log.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/mesh/mesh_continuum/grid_face_histogram.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"
#include "framework/utils/timer.h"
#include <cmath>
#include <stdexcept>

namespace opensn
{

#ifndef __OPENSN_WITH_GPU__
namespace detail
{

void
BuildAAHGPUFludsCommonData(SweepRuntime& runtime,
                           const SpatialDiscretization& discretization,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  static_cast<void>(runtime);
  static_cast<void>(discretization);
  static_cast<void>(grid_nodal_mappings);
  throw std::runtime_error("BuildAAHGPUFludsCommonData: OPENSN_WITH_CUDA not enabled.");
}

void
BuildCBCGPUFludsCommonData(SweepRuntime& runtime,
                           const SpatialDiscretization& discretization,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  static_cast<void>(runtime);
  static_cast<void>(discretization);
  static_cast<void>(grid_nodal_mappings);
  throw std::runtime_error("BuildCBCGPUFludsCommonData: OPENSN_WITH_CUDA not enabled.");
}

} // namespace detail
#endif

namespace
{

// Developer-only diagnostic for tiny sweep graphs. Populate with AAH_SPDS ids when debugging.
const std::vector<int> SWEEP_ORDER_DIRECTIONS_TO_PRINT = {};

void
AppendNonEmptyGrouping(UniqueSOGroupings& unique_so_groupings, DirIDs dir_ids)
{
  if (not dir_ids.empty())
    unique_so_groupings.push_back(std::move(dir_ids));
}

DirIDToSOMap
BuildDirectionToSweepOrderingMap(const UniqueSOGroupings& unique_so_groupings)
{
  DirIDToSOMap dir_id_to_so_map;
  for (size_t so_grouping_id = 0; so_grouping_id < unique_so_groupings.size(); ++so_grouping_id)
    for (const size_t dir_id : unique_so_groupings[so_grouping_id])
      dir_id_to_so_map[dir_id] = so_grouping_id;

  return dir_id_to_so_map;
}

std::pair<UniqueSOGroupings, DirIDToSOMap>
AssociateSOsAndDirections(const std::string& problem_name,
                          const std::shared_ptr<MeshContinuum>& grid,
                          const AngularQuadrature& quadrature,
                          AngleAggregationType agg_type,
                          GeometryType geometry_type)
{
  if (quadrature.GetOmegas().empty())
    throw std::logic_error(problem_name + ": Quadrature with no omegas cannot be used");
  if (quadrature.GetWeights().empty())
    throw std::logic_error(problem_name + ": Quadrature with no weights cannot be used");

  UniqueSOGroupings unique_so_groupings;
  switch (agg_type)
  {
    case AngleAggregationType::SINGLE:
    {
      if (geometry_type == GeometryType::TWOD_CYLINDRICAL)
      {
        const auto* product_quad = dynamic_cast<const ProductQuadrature*>(&quadrature);
        if (product_quad)
        {
          for (const auto& dir_set : product_quad->GetDirectionMap())
            for (const auto dir_id : dir_set.second)
              AppendNonEmptyGrouping(unique_so_groupings, {dir_id});
        }
        else
        {
          const size_t num_dirs = quadrature.GetNumAngles();
          for (size_t n = 0; n < num_dirs; ++n)
            AppendNonEmptyGrouping(unique_so_groupings, {n});
        }
      }
      else
      {
        const size_t num_dirs = quadrature.GetNumAngles();
        for (size_t n = 0; n < num_dirs; ++n)
          AppendNonEmptyGrouping(unique_so_groupings, {n});
      }
      break;
    }
    case AngleAggregationType::POLAR:
    {
      if (grid->GetType() != ORTHOGONAL and grid->GetDimension() != 2 and not grid->Extruded())
        throw std::logic_error(
          problem_name +
          ": The simulation is using polar angle aggregation for which only certain geometry "
          "types are supported, i.e., ORTHOGONAL, 2D or 3D EXTRUDED");

      const auto quad_type = quadrature.GetType();
      if (quad_type != AngularQuadratureType::PRODUCT_QUADRATURE)
        throw std::logic_error(problem_name +
                               ": The simulation is using polar angle aggregation for which only "
                               "Product-type quadratures are supported");

      try
      {
        const auto& product_quad = dynamic_cast<const ProductQuadrature&>(quadrature);

        const auto& azimuthal_angles = product_quad.GetAzimuthalAngles();
        const auto& polar_angles = product_quad.GetPolarAngles();
        const auto num_azimuthal = azimuthal_angles.size();
        const auto num_polar = polar_angles.size();

        std::vector<size_t> upward_polar_ids;
        std::vector<size_t> downward_polar_ids;
        for (size_t p = 0; p < num_polar; ++p)
          if (polar_angles[p] > M_PI_2)
            upward_polar_ids.push_back(p);
          else
            downward_polar_ids.push_back(p);

        auto MapPolarAndAzimuthalIDs =
          [&product_quad, &unique_so_groupings](const DirIDs& polar_ids, size_t azimuthal_id)
        {
          DirIDs dir_ids;
          dir_ids.reserve(polar_ids.size());
          for (const size_t p : polar_ids)
            dir_ids.push_back(product_quad.GetAngleNum(p, azimuthal_id));
          AppendNonEmptyGrouping(unique_so_groupings, std::move(dir_ids));
        };

        for (size_t a = 0; a < num_azimuthal; ++a)
        {
          if (not upward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(upward_polar_ids, a);
          if (not downward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(downward_polar_ids, a);
        }
      }
      catch (const std::bad_cast&)
      {
        throw std::runtime_error(problem_name +
                                 ": Casting the angular quadrature to the product quadrature base "
                                 "failed");
      }

      break;
    }
    case AngleAggregationType::AZIMUTHAL:
    {
      if (geometry_type != GeometryType::ONED_SPHERICAL and
          geometry_type != GeometryType::TWOD_CYLINDRICAL)
        throw std::logic_error(problem_name +
                               ": AZIMUTHAL aggregation is only valid for TWOD_CYLINDRICAL "
                               "geometry");

      const auto quad_type = quadrature.GetType();
      if (quad_type != AngularQuadratureType::PRODUCT_QUADRATURE)
        throw std::logic_error(problem_name +
                               ": AZIMUTHAL aggregation is only valid for TWOD_CYLINDRICAL "
                               "geometry.");

      try
      {
        const auto& product_quad = dynamic_cast<const ProductQuadrature&>(quadrature);

        for (const auto& dir_set : product_quad.GetDirectionMap())
        {
          std::vector<unsigned int> group1;
          std::vector<unsigned int> group2;
          for (const auto& dir_id : dir_set.second)
            if (quadrature.GetAbscissa(dir_id).phi > M_PI_2)
              group1.push_back(dir_id);
            else
              group2.push_back(dir_id);

          AppendNonEmptyGrouping(unique_so_groupings, {group1.begin(), group1.end()});
          AppendNonEmptyGrouping(unique_so_groupings, {group2.begin(), group2.end()});
        }
      }
      catch (const std::bad_cast&)
      {
        throw std::runtime_error(problem_name +
                                 ": Casting the angular quadrature to the product quadrature base "
                                 "failed");
      }

      break;
    }
    default:
      throw std::invalid_argument(problem_name + ": Called with UNDEFINED angle aggregation type");
  }

  return {unique_so_groupings, BuildDirectionToSweepOrderingMap(unique_so_groupings)};
}

void
BuildSweepOrderingGroups(SweepRuntime& runtime,
                         const std::string& problem_name,
                         const std::vector<LBSGroupset>& groupsets,
                         const std::shared_ptr<MeshContinuum>& grid,
                         GeometryType geometry_type,
                         std::map<std::shared_ptr<AngularQuadrature>, bool>& allow_cycles_map)
{
  for (const auto& groupset : groupsets)
  {
    if (runtime.quadrature_unq_so_grouping_map.count(groupset.quadrature) == 0)
    {
      runtime.quadrature_unq_so_grouping_map[groupset.quadrature] = AssociateSOsAndDirections(
        problem_name, grid, *groupset.quadrature, groupset.angleagg_method, geometry_type);
    }

    if (allow_cycles_map.count(groupset.quadrature) == 0)
      allow_cycles_map[groupset.quadrature] = groupset.allow_cycles;
  }
}

template <typename CreateSPDS>
void
BuildSPDSForSweepOrderings(SweepRuntime& runtime, CreateSPDS create_spds)
{
  for (const auto& [quadrature, info] : runtime.quadrature_unq_so_grouping_map)
  {
    auto& spds_list = runtime.quadrature_spds_map[quadrature];
    const auto& unique_so_groupings = info.first;
    for (const auto& so_grouping : unique_so_groupings)
    {
      const size_t master_dir_id = so_grouping.front();
      const auto& omega = quadrature->GetOmega(master_dir_id);
      spds_list.push_back(create_spds(quadrature, spds_list.size(), omega));
    }
  }
}

void
BuildAAHSPDS(SweepRuntime& runtime,
             const std::shared_ptr<MeshContinuum>& grid,
             const std::map<std::shared_ptr<AngularQuadrature>, bool>& allow_cycles_map,
             bool use_gpus)
{
  log.Log0Verbose1() << program_timer.GetTimeString() << " Initializing AAH SPDS.";
  BuildSPDSForSweepOrderings(
    runtime,
    [&grid, &allow_cycles_map, use_gpus](const std::shared_ptr<AngularQuadrature>& quadrature,
                                         size_t sweep_ordering_id,
                                         const Vector3& omega) -> std::shared_ptr<SPDS>
    {
      return std::make_shared<AAH_SPDS>(static_cast<int>(sweep_ordering_id),
                                        omega,
                                        grid,
                                        allow_cycles_map.at(quadrature),
                                        use_gpus);
    });
}

void
BuildCBCSPDS(SweepRuntime& runtime,
             const std::shared_ptr<MeshContinuum>& grid,
             const std::map<std::shared_ptr<AngularQuadrature>, bool>& allow_cycles_map)
{
  BuildSPDSForSweepOrderings(
    runtime,
    [&grid, &allow_cycles_map](const std::shared_ptr<AngularQuadrature>& quadrature,
                               size_t,
                               const Vector3& omega) -> std::shared_ptr<SPDS>
    { return std::make_shared<CBC_SPDS>(omega, grid, allow_cycles_map.at(quadrature)); });
}

std::vector<std::shared_ptr<AAH_SPDS>>
GetAAHSPDSList(const SweepRuntime& runtime)
{
  std::vector<std::shared_ptr<AAH_SPDS>> aah_spds_list;
  // The flattened order is the MPI serialization key for gathered AAH sweep FAS records. Keep this
  // traversal deterministic across ranks and in sync with ApplyAAHSweepFAS.
  for (const auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    for (const auto& spds : spds_list)
      aah_spds_list.push_back(std::static_pointer_cast<AAH_SPDS>(spds));

  return aah_spds_list;
}

int
GetSweepGraphOwner(size_t spds_ordinal)
{
  return static_cast<int>(spds_ordinal % static_cast<size_t>(opensn::mpi_comm.size()));
}

void
AccumulateAAHGlobalEdgeWeights(const std::vector<std::shared_ptr<AAH_SPDS>>& spds_list)
{
  const int comm_size = opensn::mpi_comm.size();
  const int matrix_size = comm_size * comm_size;
  std::vector<int> recv_counts(comm_size, comm_size);
  std::vector<int> recv_displacements(comm_size, 0);
  for (int loc = 0; loc < comm_size; ++loc)
    recv_displacements[loc] = loc * comm_size;

  for (size_t spds_ordinal = 0; spds_ordinal < spds_list.size(); ++spds_ordinal)
  {
    const auto& spds = spds_list[spds_ordinal];
    const int owner = GetSweepGraphOwner(spds_ordinal);

    const auto local_row = spds->ComputeLocalLocationEdgeWeights();
    std::vector<double> recv;
    if (opensn::mpi_comm.rank() == owner)
      recv.assign(matrix_size, 0.0);
    opensn::mpi_comm.gather(local_row, recv, recv_counts, recv_displacements, owner);

    if (opensn::mpi_comm.rank() == owner)
      spds->SetGlobalEdgeWeights(std::move(recv));
  }
}

void
BuildOwnedAAHSweepFAS(const std::vector<std::shared_ptr<AAH_SPDS>>& spds_list)
{
  log.Log0Verbose1() << program_timer.GetTimeString() << " Build global sweep FAS for each SPDS.";
  for (size_t spds_ordinal = 0; spds_ordinal < spds_list.size(); ++spds_ordinal)
    if (opensn::mpi_comm.rank() == GetSweepGraphOwner(spds_ordinal))
      spds_list[spds_ordinal]->BuildGlobalSweepFAS();
}

std::vector<int>
GatherAAHSweepFAS(const std::vector<std::shared_ptr<AAH_SPDS>>& spds_list)
{
  log.Log0Verbose1() << program_timer.GetTimeString() << " Gather FAS for each SPDS.";
  std::vector<int> local_edges_to_remove;
  for (size_t spds_ordinal = 0; spds_ordinal < spds_list.size(); ++spds_ordinal)
  {
    if (opensn::mpi_comm.rank() == GetSweepGraphOwner(spds_ordinal))
    {
      const auto& edges_to_remove = spds_list[spds_ordinal]->GetGlobalSweepFAS();
      local_edges_to_remove.push_back(static_cast<int>(spds_ordinal));
      local_edges_to_remove.push_back(static_cast<int>(edges_to_remove.size()));
      local_edges_to_remove.insert(
        local_edges_to_remove.end(), edges_to_remove.begin(), edges_to_remove.end());
    }
  }

  int local_size = static_cast<int>(local_edges_to_remove.size());
  std::vector<int> receive_counts(opensn::mpi_comm.size(), 0);
  std::vector<int> displacements(opensn::mpi_comm.size(), 0);
  mpi_comm.all_gather(local_size, receive_counts);

  int total_size = 0;
  for (size_t i = 0; i < receive_counts.size(); ++i)
  {
    displacements[i] = total_size;
    total_size += receive_counts[i];
  }

  std::vector<int> global_edges_to_remove(total_size, 0);
  mpi_comm.all_gather(local_edges_to_remove, global_edges_to_remove, receive_counts, displacements);

  return global_edges_to_remove;
}

void
ApplyAAHSweepFAS(const std::vector<std::shared_ptr<AAH_SPDS>>& spds_list,
                 const std::vector<int>& global_edges_to_remove)
{
  size_t offset = 0;
  while (offset < global_edges_to_remove.size())
  {
    if (offset + 2 > global_edges_to_remove.size())
      OpenSnLogicalError("Malformed AAH sweep FAS payload.");

    const auto spds_ordinal = static_cast<size_t>(global_edges_to_remove[offset++]);
    const auto num_edges = global_edges_to_remove[offset++];
    if (num_edges < 0 or offset + static_cast<size_t>(num_edges) > global_edges_to_remove.size())
      OpenSnLogicalError("Malformed AAH sweep FAS payload.");

    std::vector<int> edges;
    edges.reserve(num_edges);
    for (int i = 0; i < num_edges; ++i)
      edges.emplace_back(global_edges_to_remove[offset++]);

    if (spds_ordinal >= spds_list.size())
      OpenSnLogicalError("Invalid SPDS ordinal gathered while building AAH sweep graph.");
    spds_list[spds_ordinal]->SetGlobalSweepFAS(std::move(edges));
  }
}

void
BuildAAHSweepTDGs(const std::vector<std::shared_ptr<AAH_SPDS>>& spds_list)
{
  log.Log0Verbose1() << program_timer.GetTimeString() << " Build global sweep TDGs.";
  for (const auto& spds : spds_list)
    spds->BuildGlobalSweepTDG();
}

void
PrintRequestedSweepGraphs(const std::vector<std::shared_ptr<AAH_SPDS>>& spds_list)
{
  for (const auto& spds : spds_list)
    for (const int dir_id : SWEEP_ORDER_DIRECTIONS_TO_PRINT)
      if (spds->GetId() == dir_id)
        spds->PrintGhostedGraph();
}

void
BuildAAHGlobalSweepGraph(SweepRuntime& runtime)
{
  // Creating an AAH SPDS can be expensive, so the global graph work is split across MPI ranks:
  // accumulate graph edge weights, build the feedback arc set on owning ranks, gather the FAS on
  // all ranks, then build each global sweep task dependency graph locally.
  const auto spds_list = GetAAHSPDSList(runtime);
  AccumulateAAHGlobalEdgeWeights(spds_list);
  BuildOwnedAAHSweepFAS(spds_list);
  ApplyAAHSweepFAS(spds_list, GatherAAHSweepFAS(spds_list));
  BuildAAHSweepTDGs(spds_list);
  PrintRequestedSweepGraphs(spds_list);
}

void
BuildAAHCPUFludsCommonData(SweepRuntime& runtime,
                           const std::shared_ptr<MeshContinuum>& grid,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  const auto grid_face_histogram = grid->MakeGridFaceHistogram();
  for (const auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    for (const auto& spds : spds_list)
      runtime.quadrature_fluds_commondata_map[quadrature].push_back(
        std::make_unique<AAH_FLUDSCommonData>(grid_nodal_mappings, *spds, *grid_face_histogram));
}

void
BuildCBCCPUFludsCommonData(SweepRuntime& runtime,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  for (const auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    for (const auto& spds : spds_list)
      runtime.quadrature_fluds_commondata_map[quadrature].push_back(
        std::make_unique<CBC_FLUDSCommonData>(*spds, grid_nodal_mappings));
}

} // namespace

SweepRuntime
BuildSweepRuntime(const std::string& problem_name,
                  const std::vector<LBSGroupset>& groupsets,
                  const std::shared_ptr<MeshContinuum>& grid,
                  const std::string& sweep_type,
                  bool use_gpus,
                  const SpatialDiscretization& discretization,
                  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  SweepRuntime runtime;
  std::map<std::shared_ptr<AngularQuadrature>, bool> quadrature_allow_cycles_map;

  const auto geometry_type = grid->GetGeometryType();
  BuildSweepOrderingGroups(
    runtime, problem_name, groupsets, grid, geometry_type, quadrature_allow_cycles_map);

  if (sweep_type == "AAH")
  {
    BuildAAHSPDS(runtime, grid, quadrature_allow_cycles_map, use_gpus);
    BuildAAHGlobalSweepGraph(runtime);
  }
  else if (sweep_type == "CBC")
    BuildCBCSPDS(runtime, grid, quadrature_allow_cycles_map);
  else
    OpenSnInvalidArgument("Unsupported sweep type \"" + sweep_type + "\"");

  opensn::mpi_comm.barrier();

  if (not use_gpus)
  {
    if (sweep_type == "AAH")
      BuildAAHCPUFludsCommonData(runtime, grid, grid_nodal_mappings);
    else if (sweep_type == "CBC")
      BuildCBCCPUFludsCommonData(runtime, grid_nodal_mappings);
  }
  else
  {
    if (sweep_type == "AAH")
      detail::BuildAAHGPUFludsCommonData(runtime, discretization, grid_nodal_mappings);
    else if (sweep_type == "CBC")
      detail::BuildCBCGPUFludsCommonData(runtime, discretization, grid_nodal_mappings);
  }

  return runtime;
}

} // namespace opensn
