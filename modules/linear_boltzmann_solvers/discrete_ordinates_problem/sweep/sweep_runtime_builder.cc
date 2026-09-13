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
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep_parallel_for.h"
#include "framework/utils/error.h"
#include "framework/utils/timer.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>

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

template <typename Function>
void
RunAAHSetupCollectively(const Function& function)
{
  bool local_failed = false;
  std::array<char, 2048> error_buffer{};
  auto CaptureError = [&](std::string_view message)
  {
    local_failed = true;
    if (message.empty())
      message = "Unknown exception during AAH setup.";
    std::copy_n(
      message.begin(), std::min(message.size(), error_buffer.size() - 1), error_buffer.begin());
  };

  try
  {
    function();
  }
  catch (const std::exception& error)
  {
    CaptureError(error.what());
  }
  catch (...)
  {
    CaptureError("Unknown exception during AAH setup.");
  }

  const int comm_size = opensn::mpi_comm.size();
  const int local_failure_rank = local_failed ? opensn::mpi_comm.rank() : comm_size;
  int failure_rank = comm_size;
  opensn::mpi_comm.all_reduce(local_failure_rank, failure_rank, mpi::op::min<int>());
  if (failure_rank == comm_size)
    return;

  opensn::mpi_comm.broadcast(
    error_buffer.data(), static_cast<int>(error_buffer.size()), failure_rank);
  throw std::logic_error("AAH sweep setup failed on rank " + std::to_string(failure_rank) + ": " +
                         error_buffer.data());
}

template <typename Function>
void
RunAAHCollectiveSequence(const Function& function)
{
  try
  {
    function();
  }
  catch (const std::exception& error)
  {
    log.LogAllError() << "AAH setup failed during MPI exchange: " << error.what();
    opensn::mpi_comm.abort(1);
    throw;
  }
  catch (...)
  {
    opensn::mpi_comm.abort(1);
    throw;
  }
}

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
BuildSweepOrderingGroups(
  SweepRuntime& runtime,
  const std::string& problem_name,
  const std::vector<LBSGroupset>& groupsets,
  const std::shared_ptr<MeshContinuum>& grid,
  GeometryType geometry_type,
  std::map<std::shared_ptr<AngularQuadrature>, bool>& allow_cycles_map,
  std::map<std::shared_ptr<AngularQuadrature>, AngleAggregationType>& angle_aggregation_map)
{
  for (const auto& groupset : groupsets)
  {
    if (runtime.quadrature_unq_so_grouping_map.count(groupset.quadrature) == 0)
    {
      runtime.quadrature_order.push_back(groupset.quadrature);
      runtime.quadrature_unq_so_grouping_map[groupset.quadrature] = AssociateSOsAndDirections(
        problem_name, grid, *groupset.quadrature, groupset.angleagg_method, geometry_type);
    }

    if (allow_cycles_map.count(groupset.quadrature) == 0)
      allow_cycles_map[groupset.quadrature] = groupset.allow_cycles;
    else if (allow_cycles_map.at(groupset.quadrature) != groupset.allow_cycles)
      throw std::invalid_argument(problem_name +
                                  ": Groupsets sharing a quadrature must use the same "
                                  "allow_cycles setting.");

    if (angle_aggregation_map.count(groupset.quadrature) == 0)
      angle_aggregation_map[groupset.quadrature] = groupset.angleagg_method;
    else if (angle_aggregation_map.at(groupset.quadrature) != groupset.angleagg_method)
      throw std::invalid_argument(problem_name +
                                  ": Groupsets sharing a quadrature must use the same angle "
                                  "aggregation type.");
  }
}

template <typename CreateSPDS>
void
BuildSPDSForSweepOrderings(SweepRuntime& runtime, CreateSPDS create_spds)
{
  struct WorkItem
  {
    std::shared_ptr<AngularQuadrature> quadrature;
    size_t sweep_ordering_id = 0;
    Vector3 omega;
  };
  std::vector<WorkItem> work;
  for (const auto& quadrature : runtime.quadrature_order)
  {
    const auto& info = runtime.quadrature_unq_so_grouping_map.at(quadrature);
    const auto& unique_so_groupings = info.first;
    size_t sweep_ordering_id = 0;
    for (const auto& so_grouping : unique_so_groupings)
      work.push_back({quadrature, sweep_ordering_id++, quadrature->GetOmega(so_grouping.front())});
  }

  for (auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    spds_list.clear();
  if (work.empty())
    return;
  std::vector<std::shared_ptr<SPDS>> result(work.size());
  const size_t nthreads = std::min<size_t>(std::max(1U, opensn_num_threads), work.size());
  log.Log() << program_timer.GetTimeString() << " SPDS construction: " << work.size() << " angles ("
            << nthreads << " threads(s)).";
  ParallelFor(
    work.size(),
    nthreads,
    [&](size_t i)
    { result[i] = create_spds(work[i].quadrature, work[i].sweep_ordering_id, work[i].omega); });
  for (size_t i = 0; i < work.size(); ++i)
    runtime.quadrature_spds_map[work[i].quadrature].push_back(std::move(result[i]));
  log.Log() << program_timer.GetTimeString() << " SPDS construction done.";
}

std::vector<std::shared_ptr<AAH_SPDS>> GetAAHSPDSList(const SweepRuntime& runtime);

void
BuildAAHSPDS(SweepRuntime& runtime,
             const std::shared_ptr<MeshContinuum>& grid,
             const SPDSFaceNeighborInfoVec& face_neighbor_info,
             const std::map<std::shared_ptr<AngularQuadrature>, bool>& allow_cycles_map,
             bool use_gpus)
{
  BuildSPDSForSweepOrderings(runtime,
                             [&grid, &face_neighbor_info, &allow_cycles_map](
                               const std::shared_ptr<AngularQuadrature>& quadrature,
                               size_t sweep_ordering_id,
                               const Vector3& omega) -> std::shared_ptr<SPDS>
                             {
                               return std::make_shared<AAH_SPDS>(
                                 static_cast<int>(sweep_ordering_id),
                                 omega,
                                 grid,
                                 face_neighbor_info,
                                 allow_cycles_map.at(quadrature),
                                 false);
                             });

  const auto spds_list = GetAAHSPDSList(runtime);
  if (use_gpus)
    for (const auto& spds : spds_list)
      spds->CopySPLSDataOnDevice();
}

void
BuildCBCSPDS(SweepRuntime& runtime,
             const std::shared_ptr<MeshContinuum>& grid,
             const SPDSFaceNeighborInfoVec& face_neighbor_info,
             const std::map<std::shared_ptr<AngularQuadrature>, bool>& allow_cycles_map)
{
  BuildSPDSForSweepOrderings(runtime,
                             [&grid, &face_neighbor_info, &allow_cycles_map](
                               const std::shared_ptr<AngularQuadrature>& quadrature,
                               size_t,
                               const Vector3& omega) -> std::shared_ptr<SPDS>
                             {
                               return std::make_shared<CBC_SPDS>(
                                 omega, grid, face_neighbor_info, allow_cycles_map.at(quadrature));
                             });
}

std::vector<std::shared_ptr<AAH_SPDS>>
GetAAHSPDSList(const SweepRuntime& runtime)
{
  std::vector<std::shared_ptr<AAH_SPDS>> aah_spds_list;
  for (const auto& quadrature : runtime.quadrature_order)
  {
    const auto& spds_list = runtime.quadrature_spds_map.at(quadrature);
    for (const auto& spds : spds_list)
      aah_spds_list.push_back(std::static_pointer_cast<AAH_SPDS>(spds));
  }

  return aah_spds_list;
}

int
GetSweepGraphOwner(size_t spds_ordinal, size_t num_spds)
{
  OpenSnLogicalErrorIf(num_spds == 0, "Cannot assign an owner without AAH sweep graphs.");
  const auto comm_size = static_cast<size_t>(opensn::mpi_comm.size());
  if (num_spds <= comm_size)
    return static_cast<int>((spds_ordinal * comm_size) / num_spds);
  return static_cast<int>(spds_ordinal % comm_size);
}

void
GatherAAHGlobalEdgeWeights(const std::shared_ptr<AAH_SPDS>& spds, int owner)
{
  const int comm_size = opensn::mpi_comm.size();
  const bool is_owner = opensn::mpi_comm.rank() == owner;

  const auto local_edges = spds->ComputeLocalLocationEdgeWeights();
  OpenSnLogicalErrorIf(local_edges.size() > static_cast<size_t>(std::numeric_limits<int>::max()),
                       "Too many local AAH sweep graph edges.");
  std::vector<int> local_to_locs(local_edges.size());
  std::vector<double> local_weights(local_edges.size());
  for (size_t i = 0; i < local_edges.size(); ++i)
  {
    local_to_locs[i] = local_edges[i].first;
    local_weights[i] = local_edges[i].second;
  }

  std::vector<int> counts(comm_size, 0);
  opensn::mpi_comm.gather(static_cast<int>(local_edges.size()), counts, owner);

  std::vector<int> offsets(comm_size, 0);
  if (is_owner)
  {
    for (int r = 1; r < comm_size; ++r)
    {
      OpenSnLogicalErrorIf(counts[r - 1] > std::numeric_limits<int>::max() - offsets[r - 1],
                           "AAH sweep graph edge count exceeds the MPI count limit.");
      offsets[r] = offsets[r - 1] + counts[r - 1];
    }
    OpenSnLogicalErrorIf(counts.back() > std::numeric_limits<int>::max() - offsets.back(),
                         "AAH sweep graph edge count exceeds the MPI count limit.");
  }

  std::vector<int> recv_to_locs;
  std::vector<double> recv_weights;
  opensn::mpi_comm.gather(local_to_locs, recv_to_locs, counts, offsets, owner);
  opensn::mpi_comm.gather(local_weights, recv_weights, counts, offsets, owner);

  if (is_owner)
  {
    std::vector<AAH_SPDS::GlobalSweepEdge> global_edges;
    global_edges.reserve(recv_to_locs.size());
    for (int r = 0; r < comm_size; ++r)
      for (int i = offsets[r]; i < offsets[r] + counts[r]; ++i)
        global_edges.push_back({r, recv_to_locs[i], recv_weights[i]});
    spds->SetGlobalEdges(std::move(global_edges));
  }
}

void
DistributeAAHGlobalSweepMetadata(
  const std::vector<std::shared_ptr<AAH_SPDS>>& spds_list,
  size_t batch_begin,
  size_t batch_end,
  const std::vector<std::optional<AAH_SPDS::GlobalSweepMetadata>>& metadata)
{
  const int comm_size = opensn::mpi_comm.size();
  std::vector<int> send_counts(comm_size, 0);
  std::vector<int> send_offsets(comm_size, 0);
  std::vector<int> send_buffer;

  for (size_t i = batch_begin; i < batch_end; ++i)
  {
    const auto& entry = metadata[i - batch_begin];
    if (not entry.has_value())
      continue;
    const auto& graph_metadata = *entry;
    OpenSnLogicalErrorIf(
      graph_metadata.location_depths.size() != static_cast<size_t>(comm_size) or
        graph_metadata.delayed_dependencies.size() != static_cast<size_t>(comm_size) or
        graph_metadata.delayed_successors.size() != static_cast<size_t>(comm_size),
      "Invalid AAH sweep graph metadata size.");

    for (int r = 0; r < comm_size; ++r)
    {
      const size_t record_size = 4 + graph_metadata.delayed_dependencies[r].size() +
                                 graph_metadata.delayed_successors[r].size();
      OpenSnLogicalErrorIf(record_size > static_cast<size_t>(std::numeric_limits<int>::max()) or
                             static_cast<size_t>(send_counts[r]) >
                               static_cast<size_t>(std::numeric_limits<int>::max()) - record_size,
                           "AAH sweep metadata exceeds the MPI count limit.");
      send_counts[r] += static_cast<int>(record_size);
    }
  }

  for (int r = 1; r < comm_size; ++r)
  {
    OpenSnLogicalErrorIf(send_counts[r - 1] > std::numeric_limits<int>::max() - send_offsets[r - 1],
                         "AAH sweep metadata exceeds the MPI count limit.");
    send_offsets[r] = send_offsets[r - 1] + send_counts[r - 1];
  }
  OpenSnLogicalErrorIf(send_counts.back() > std::numeric_limits<int>::max() - send_offsets.back(),
                       "AAH sweep metadata exceeds the MPI count limit.");
  send_buffer.reserve(static_cast<size_t>(send_offsets.back()) + send_counts.back());

  for (int r = 0; r < comm_size; ++r)
    for (size_t i = batch_begin; i < batch_end; ++i)
    {
      const auto& entry = metadata[i - batch_begin];
      if (not entry.has_value())
        continue;
      const auto& graph_metadata = *entry;
      const auto& dependencies = graph_metadata.delayed_dependencies[r];
      const auto& successors = graph_metadata.delayed_successors[r];
      send_buffer.push_back(static_cast<int>(i - batch_begin));
      send_buffer.push_back(graph_metadata.location_depths[r]);
      send_buffer.push_back(static_cast<int>(dependencies.size()));
      send_buffer.push_back(static_cast<int>(successors.size()));
      send_buffer.insert(send_buffer.end(), dependencies.begin(), dependencies.end());
      send_buffer.insert(send_buffer.end(), successors.begin(), successors.end());
    }

  std::vector<int> receive_counts(comm_size, 0);
  MPI_CHECK(MPI_Alltoall(send_counts.data(),
                         1,
                         MPI_INT,
                         receive_counts.data(),
                         1,
                         MPI_INT,
                         static_cast<MPI_Comm>(mpi_comm)));

  std::vector<int> receive_offsets(comm_size, 0);
  for (int r = 1; r < comm_size; ++r)
  {
    OpenSnLogicalErrorIf(receive_counts[r - 1] >
                           std::numeric_limits<int>::max() - receive_offsets[r - 1],
                         "AAH sweep metadata exceeds the MPI count limit.");
    receive_offsets[r] = receive_offsets[r - 1] + receive_counts[r - 1];
  }
  OpenSnLogicalErrorIf(receive_counts.back() >
                         std::numeric_limits<int>::max() - receive_offsets.back(),
                       "AAH sweep metadata exceeds the MPI count limit.");
  std::vector<int> receive_buffer(static_cast<size_t>(receive_offsets.back()) +
                                  receive_counts.back());
  MPI_CHECK(MPI_Alltoallv(send_buffer.data(),
                          send_counts.data(),
                          send_offsets.data(),
                          MPI_INT,
                          receive_buffer.data(),
                          receive_counts.data(),
                          receive_offsets.data(),
                          MPI_INT,
                          static_cast<MPI_Comm>(mpi_comm)));

  std::vector<bool> received(batch_end - batch_begin, false);
  for (int source = 0; source < comm_size; ++source)
  {
    size_t offset = receive_offsets[source];
    const size_t end = offset + receive_counts[source];
    while (offset < end)
    {
      OpenSnLogicalErrorIf(end - offset < 4, "Malformed AAH sweep metadata record.");
      const int batch_offset = receive_buffer[offset++];
      const int location_depth = receive_buffer[offset++];
      const int num_dependencies = receive_buffer[offset++];
      const int num_successors = receive_buffer[offset++];
      OpenSnLogicalErrorIf(batch_offset < 0 or
                             static_cast<size_t>(batch_offset) >= received.size() or
                             received[batch_offset] or num_dependencies < 0 or num_successors < 0 or
                             static_cast<size_t>(num_dependencies) > end - offset or
                             static_cast<size_t>(num_successors) >
                               end - offset - static_cast<size_t>(num_dependencies),
                           "Malformed AAH sweep metadata record.");

      const auto record_begin =
        receive_buffer.begin() + static_cast<std::vector<int>::difference_type>(offset);
      const auto dependency_end = record_begin + num_dependencies;
      const auto successor_end = dependency_end + num_successors;
      spds_list[batch_begin + batch_offset]->SetGlobalSweepMetadata(
        location_depth, {record_begin, dependency_end}, {dependency_end, successor_end});
      received[batch_offset] = true;
      offset += static_cast<size_t>(num_dependencies) + num_successors;
    }
  }
  OpenSnLogicalErrorIf(std::find(received.begin(), received.end(), false) != received.end(),
                       "Missing AAH sweep metadata record.");
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
  const auto spds_list = GetAAHSPDSList(runtime);
  const auto batch_size = static_cast<size_t>(opensn::mpi_comm.size());
  for (size_t batch_begin = 0; batch_begin < spds_list.size(); batch_begin += batch_size)
  {
    const size_t batch_end = std::min(batch_begin + batch_size, spds_list.size());
    for (size_t i = batch_begin; i < batch_end; ++i)
      RunAAHCollectiveSequence(
        [&] { GatherAAHGlobalEdgeWeights(spds_list[i], GetSweepGraphOwner(i, spds_list.size())); });

    std::vector<std::optional<AAH_SPDS::GlobalSweepMetadata>> metadata;
    RunAAHSetupCollectively(
      [&]
      {
        metadata.resize(batch_end - batch_begin);
        for (size_t i = batch_begin; i < batch_end; ++i)
          if (opensn::mpi_comm.rank() == GetSweepGraphOwner(i, spds_list.size()))
            metadata[i - batch_begin] = spds_list[i]->BuildGlobalSweepMetadata();
      });

    RunAAHCollectiveSequence(
      [&] { DistributeAAHGlobalSweepMetadata(spds_list, batch_begin, batch_end, metadata); });
  }
  PrintRequestedSweepGraphs(spds_list);
}

void
BuildAAHCPUFludsCommonData(SweepRuntime& runtime,
                           const std::shared_ptr<MeshContinuum>& grid,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  struct WorkItem
  {
    std::shared_ptr<AngularQuadrature> quadrature;
    std::shared_ptr<SPDS> spds;
  };
  std::vector<WorkItem> work;
  for (const auto& quadrature : runtime.quadrature_order)
  {
    const auto& spds_list = runtime.quadrature_spds_map.at(quadrature);
    for (const auto& spds : spds_list)
      work.push_back({quadrature, spds});
  }
  if (work.empty())
    return;

  const auto grid_face_histogram = grid->MakeGridFaceHistogram();
  const size_t nthreads = std::min<size_t>(std::max(1U, opensn_num_threads), work.size());
  log.Log() << program_timer.GetTimeString() << " FLUDS construction: " << work.size()
            << " angles (" << nthreads << " threads(s)).";
  std::vector<std::unique_ptr<AAH_FLUDSCommonData>> result(work.size());
  ParallelFor(work.size(),
              nthreads,
              [&](size_t i)
              {
                result[i] = AAH_FLUDSCommonData::MakeAlpha(
                  grid_nodal_mappings, *work[i].spds, *grid_face_histogram);
              });

  // Inter-rank face exchange must retain identical ordering on every rank.
  for (size_t i = 0; i < work.size(); ++i)
    result[i]->FinalizeBeta(*work[i].spds);
  for (size_t i = 0; i < work.size(); ++i)
    runtime.quadrature_fluds_commondata_map[work[i].quadrature].push_back(std::move(result[i]));
  log.Log() << program_timer.GetTimeString() << " FLUDS construction done.";
}

void
BuildCBCCPUFludsCommonData(SweepRuntime& runtime,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  struct WorkItem
  {
    std::shared_ptr<AngularQuadrature> quadrature;
    std::shared_ptr<SPDS> spds;
  };
  std::vector<WorkItem> work;
  for (const auto& quadrature : runtime.quadrature_order)
  {
    const auto& spds_list = runtime.quadrature_spds_map.at(quadrature);
    for (const auto& spds : spds_list)
      work.push_back({quadrature, spds});
  }
  if (work.empty())
    return;

  const size_t nthreads = std::min<size_t>(std::max(1U, opensn_num_threads), work.size());
  log.Log() << program_timer.GetTimeString() << " FLUDS construction: " << work.size()
            << " angles (" << nthreads << " threads(s)).";
  std::vector<std::unique_ptr<CBC_FLUDSCommonData>> result(work.size());
  ParallelFor(work.size(),
              nthreads,
              [&](size_t i)
              { result[i] = CBC_FLUDSCommonData::MakeAlpha(*work[i].spds, grid_nodal_mappings); });
  for (auto& fluds : result)
    fluds->FinalizeBeta();
  for (size_t i = 0; i < work.size(); ++i)
    runtime.quadrature_fluds_commondata_map[work[i].quadrature].push_back(std::move(result[i]));
  log.Log() << program_timer.GetTimeString() << " FLUDS construction done.";
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
  std::map<std::shared_ptr<AngularQuadrature>, AngleAggregationType> quadrature_aggregation_map;

  const auto geometry_type = grid->GetGeometryType();
  BuildSweepOrderingGroups(runtime,
                           problem_name,
                           groupsets,
                           grid,
                           geometry_type,
                           quadrature_allow_cycles_map,
                           quadrature_aggregation_map);
  const auto face_neighbor_info = BuildSPDSFaceNeighborInfo(*grid);

  if (sweep_type == "AAH")
  {
    RunAAHSetupCollectively(
      [&]
      { BuildAAHSPDS(runtime, grid, face_neighbor_info, quadrature_allow_cycles_map, use_gpus); });
    BuildAAHGlobalSweepGraph(runtime);
  }
  else if (sweep_type == "CBC")
    BuildCBCSPDS(runtime, grid, face_neighbor_info, quadrature_allow_cycles_map);
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
