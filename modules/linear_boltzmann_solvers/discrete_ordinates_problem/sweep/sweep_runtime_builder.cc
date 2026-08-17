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
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>

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

// Pass any exception message to `on_failure` from within its catch block.
template <typename Function, typename OnFailure>
void
RunAAHSetupStep(const Function& function, const OnFailure& on_failure)
{
  try
  {
    function();
  }
  catch (const std::exception& error)
  {
    std::string_view message = error.what();
    on_failure(message.empty() ? std::string_view("Unknown exception during AAH setup.") : message);
  }
  catch (...)
  {
    on_failure(std::string_view("Unknown exception during AAH setup."));
  }
}

// Propagate failures collectively from work that does not call MPI.
template <typename Function>
void
RunAAHSetupCollectively(const Function& function)
{
  bool local_failed = false;
  std::array<char, 2048> error_buffer{};
  RunAAHSetupStep(function,
                  [&](std::string_view message)
                  {
                    local_failed = true;
                    std::copy_n(message.begin(),
                                std::min(message.size(), error_buffer.size() - 1),
                                error_buffer.begin());
                  });

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

// Abort if work containing collectives fails because peers may already be blocked.
template <typename Function>
void
RunAAHCollectiveSequence(const Function& function)
{
  RunAAHSetupStep(function,
                  [](std::string_view message)
                  {
                    log.LogAllError() << "AAH setup failed during MPI exchange: " << message;
                    opensn::mpi_comm.abort(1);
                    throw;
                  });
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
std::vector<std::shared_ptr<CBC_SPDS>> GetCBCSPDSList(const SweepRuntime& runtime);

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
             const std::vector<LBSGroupset>& groupsets,
             const std::shared_ptr<MeshContinuum>& grid,
             const SPDSFaceNeighborInfoVec& face_neighbor_info,
             const std::map<std::shared_ptr<AngularQuadrature>, bool>& allow_cycles_map,
             bool use_gpus)
{
  struct WorkItem
  {
    std::shared_ptr<AngularQuadrature> quadrature;
    Vector3 omega;
  };

  std::vector<std::shared_ptr<AngularQuadrature>> ordered_quadratures;
  ordered_quadratures.reserve(groupsets.size());
  for (const auto& groupset : groupsets)
    if (std::ranges::find(ordered_quadratures, groupset.quadrature) == ordered_quadratures.end())
      ordered_quadratures.push_back(groupset.quadrature);

  std::vector<WorkItem> work;
  for (const auto& quadrature : ordered_quadratures)
  {
    const auto& unique_so_groupings = runtime.quadrature_unq_so_grouping_map.at(quadrature).first;
    for (const auto& so_grouping : unique_so_groupings)
      work.push_back({quadrature, quadrature->GetOmega(so_grouping.front())});
  }

  for (auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    spds_list.clear();
  if (work.empty())
    return;

  std::vector<std::shared_ptr<SPDS>> result(work.size());
  const size_t nthreads = std::min<size_t>(std::max(1U, opensn_num_threads), work.size());
  log.Log() << program_timer.GetTimeString() << " SPDS construction: " << work.size() << " angles ("
            << nthreads << " thread(s)).";
  ParallelFor(work.size(),
              nthreads,
              [&](size_t i)
              {
                result[i] = std::make_shared<CBC_SPDS>(static_cast<int>(i),
                                                       work[i].omega,
                                                       grid,
                                                       face_neighbor_info,
                                                       allow_cycles_map.at(work[i].quadrature) and
                                                         not use_gpus);
              });
  for (size_t i = 0; i < work.size(); ++i)
    runtime.quadrature_spds_map[work[i].quadrature].push_back(std::move(result[i]));
  log.Log() << program_timer.GetTimeString() << " SPDS construction done.";
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

std::vector<std::shared_ptr<CBC_SPDS>>
GetCBCSPDSList(const SweepRuntime& runtime)
{
  std::vector<std::shared_ptr<CBC_SPDS>> cbc_spds_list;
  for (const auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    for (const auto& spds : spds_list)
      cbc_spds_list.push_back(std::static_pointer_cast<CBC_SPDS>(spds));

  std::ranges::sort(cbc_spds_list, {}, &CBC_SPDS::GetId);

  return cbc_spds_list;
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

std::vector<int>
BuildOffsets(const std::vector<int>& counts, const char* error_message)
{
  std::vector<int> offsets(counts.size(), 0);
  for (size_t i = 1; i < counts.size(); ++i)
  {
    OpenSnLogicalErrorIf(counts[i - 1] > std::numeric_limits<int>::max() - offsets[i - 1],
                         error_message);
    offsets[i] = offsets[i - 1] + counts[i - 1];
  }
  OpenSnLogicalErrorIf(counts.back() > std::numeric_limits<int>::max() - offsets.back(),
                       error_message);
  return offsets;
}

void
GatherAAHGlobalEdgeWeights(const std::vector<std::shared_ptr<AAH_SPDS>>& spds_list,
                           size_t batch_begin,
                           size_t batch_end)
{
  const int comm_size = opensn::mpi_comm.size();
  const int rank = opensn::mpi_comm.rank();
  const size_t batch_size = batch_end - batch_begin;
  const char* const error_message = "AAH sweep graph edge count exceeds the MPI count limit.";

  // The deterministic owner map also gives every rank its batch ownership count.
  std::vector<int> owner(batch_size);
  int owned_count = 0;
  for (size_t local_i = 0; local_i < batch_size; ++local_i)
  {
    owner[local_i] = GetSweepGraphOwner(batch_begin + local_i, spds_list.size());
    if (owner[local_i] == rank)
      ++owned_count;
  }

  // Compute each sparse local edge list once and use it to size the packed buffers.
  std::vector<std::vector<std::pair<int, double>>> local_edges(batch_size);
  std::vector<int> send_int_counts(comm_size, 0);
  std::vector<int> send_double_counts(comm_size, 0);
  for (size_t local_i = 0; local_i < batch_size; ++local_i)
  {
    auto& edges = local_edges[local_i];
    edges = spds_list[batch_begin + local_i]->ComputeLocalLocationEdgeWeights();
    const auto num_edges = edges.size();
    OpenSnLogicalErrorIf(num_edges > static_cast<size_t>(std::numeric_limits<int>::max()),
                         "Too many local AAH sweep graph edges.");
    auto& int_count = send_int_counts[owner[local_i]];
    OpenSnLogicalErrorIf(static_cast<size_t>(int_count) >
                           static_cast<size_t>(std::numeric_limits<int>::max()) - (2 + num_edges),
                         error_message);
    int_count += static_cast<int>(2 + num_edges);
    send_double_counts[owner[local_i]] += static_cast<int>(num_edges);
  }

  const auto send_int_offsets = BuildOffsets(send_int_counts, error_message);
  const auto send_double_offsets = BuildOffsets(send_double_counts, error_message);

  // Pack each edge list as an integer header and parallel destination/weight arrays.
  std::vector<int> send_int_buffer(static_cast<size_t>(send_int_offsets.back()) +
                                   static_cast<size_t>(send_int_counts.back()));
  std::vector<double> send_double_buffer(static_cast<size_t>(send_double_offsets.back()) +
                                         static_cast<size_t>(send_double_counts.back()));
  auto int_cursor = send_int_offsets;
  auto double_cursor = send_double_offsets;
  for (size_t local_i = 0; local_i < batch_size; ++local_i)
  {
    auto& edges = local_edges[local_i];
    auto& ic = int_cursor[owner[local_i]];
    auto& dc = double_cursor[owner[local_i]];
    send_int_buffer[ic++] = static_cast<int>(local_i);
    send_int_buffer[ic++] = static_cast<int>(edges.size());
    for (const auto& [to_loc, weight] : edges)
    {
      send_int_buffer[ic++] = to_loc;
      send_double_buffer[dc++] = weight;
    }
    std::vector<std::pair<int, double>>().swap(edges);
  }

  // Derive integer receive counts from the exchanged edge counts and known record count.
  std::vector<int> receive_double_counts(comm_size, 0);
  opensn::mpi_comm.all_to_all(send_double_counts, receive_double_counts);
  const auto receive_double_offsets = BuildOffsets(receive_double_counts, error_message);

  std::vector<int> receive_int_counts(comm_size);
  for (int source = 0; source < comm_size; ++source)
    receive_int_counts[source] = 2 * owned_count + receive_double_counts[source];
  const auto receive_int_offsets = BuildOffsets(receive_int_counts, error_message);

  std::vector<int> receive_int_buffer;
  opensn::mpi_comm.all_to_all(send_int_buffer,
                              send_int_counts,
                              send_int_offsets,
                              receive_int_buffer,
                              receive_int_counts,
                              receive_int_offsets);
  std::vector<int>().swap(send_int_buffer);

  std::vector<double> receive_double_buffer;
  opensn::mpi_comm.all_to_all(send_double_buffer,
                              send_double_counts,
                              send_double_offsets,
                              receive_double_buffer,
                              receive_double_counts,
                              receive_double_offsets);
  std::vector<double>().swap(send_double_buffer);

  // Unpack, per source rank, into per-SPDS edge lists.
  std::vector<std::vector<AAH_SPDS::GlobalSweepEdge>> global_edges(batch_size);
  for (int source = 0; source < comm_size; ++source)
  {
    auto int_offset = static_cast<size_t>(receive_int_offsets[source]);
    const size_t int_end = int_offset + static_cast<size_t>(receive_int_counts[source]);
    auto double_offset = static_cast<size_t>(receive_double_offsets[source]);
    const size_t double_end = double_offset + static_cast<size_t>(receive_double_counts[source]);
    while (int_offset < int_end)
    {
      OpenSnLogicalErrorIf(int_end - int_offset < 2, "Malformed AAH sweep graph edge record.");
      const int local_i = receive_int_buffer[int_offset++];
      const int num_edges = receive_int_buffer[int_offset++];
      OpenSnLogicalErrorIf(
        local_i < 0 or std::cmp_greater_equal(local_i, batch_size) or num_edges < 0 or
          static_cast<size_t>(num_edges) > int_end - int_offset or double_offset > double_end or
          static_cast<size_t>(num_edges) > double_end - double_offset,
        "Malformed AAH sweep graph edge record.");
      auto& edges = global_edges[static_cast<size_t>(local_i)];
      edges.reserve(edges.size() + static_cast<size_t>(num_edges));
      for (int k = 0; k < num_edges; ++k)
        edges.push_back(
          {source, receive_int_buffer[int_offset++], receive_double_buffer[double_offset++]});
    }
    OpenSnLogicalErrorIf(double_offset != double_end, "Malformed AAH sweep graph edge record.");
  }

  for (size_t local_i = 0; local_i < batch_size; ++local_i)
    if (owner[local_i] == rank)
      spds_list[batch_begin + local_i]->SetGlobalEdges(std::move(global_edges[local_i]));
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

  const auto send_offsets =
    BuildOffsets(send_counts, "AAH sweep metadata exceeds the MPI count limit.");
  send_buffer.resize(static_cast<size_t>(send_offsets.back()) +
                     static_cast<size_t>(send_counts.back()));

  // Pack only locally owned metadata, preserving one contiguous region per destination.
  auto cursor = send_offsets;
  for (size_t i = batch_begin; i < batch_end; ++i)
  {
    const auto& entry = metadata[i - batch_begin];
    if (not entry.has_value())
      continue;
    const auto& graph_metadata = *entry;
    for (int r = 0; r < comm_size; ++r)
    {
      const auto& dependencies = graph_metadata.delayed_dependencies[r];
      const auto& successors = graph_metadata.delayed_successors[r];
      auto& c = cursor[r];
      send_buffer[c++] = static_cast<int>(i - batch_begin);
      send_buffer[c++] = graph_metadata.location_depths[r];
      send_buffer[c++] = static_cast<int>(dependencies.size());
      send_buffer[c++] = static_cast<int>(successors.size());
      for (const int dep : dependencies)
        send_buffer[c++] = dep;
      for (const int succ : successors)
        send_buffer[c++] = succ;
    }
  }

  std::vector<int> receive_counts(comm_size, 0);
  opensn::mpi_comm.all_to_all(send_counts, receive_counts);

  const auto receive_offsets =
    BuildOffsets(receive_counts, "AAH sweep metadata exceeds the MPI count limit.");
  std::vector<int> receive_buffer;
  opensn::mpi_comm.all_to_all(
    send_buffer, send_counts, send_offsets, receive_buffer, receive_counts, receive_offsets);

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
  const auto local_count = spds_list.size();
  const auto size_max = std::numeric_limits<size_t>::max();
  const std::array<size_t, 2> local_bounds = {local_count, size_max - local_count};
  std::array<size_t, 2> global_bounds{};
  opensn::mpi_comm.all_reduce(local_bounds.data(),
                              static_cast<int>(local_bounds.size()),
                              global_bounds.data(),
                              mpi::op::max<size_t>());
  OpenSnLogicalErrorIf(global_bounds[0] != size_max - global_bounds[1],
                       "AAH SPDS count differs across MPI ranks.");

  const auto batch_size = static_cast<size_t>(opensn::mpi_comm.size());
  for (size_t batch_begin = 0; batch_begin < spds_list.size(); batch_begin += batch_size)
  {
    const size_t batch_end = std::min(batch_begin + batch_size, spds_list.size());
    RunAAHCollectiveSequence([&]
                             { GatherAAHGlobalEdgeWeights(spds_list, batch_begin, batch_end); });

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

std::vector<std::uint64_t>
PackCBCSPDSLocationEdgeWeights(
  const std::vector<CBC_SPDS::LocationEdgeWeight>& location_edge_weights)
{
  constexpr std::size_t packed_values_per_edge = 3;
  std::vector<std::uint64_t> packed_values;
  packed_values.reserve(location_edge_weights.size() * packed_values_per_edge);
  for (const auto& edge_weight : location_edge_weights)
  {
    packed_values.push_back(static_cast<std::uint64_t>(edge_weight.upstream_location));
    packed_values.push_back(static_cast<std::uint64_t>(edge_weight.downstream_location));
    packed_values.push_back(std::bit_cast<std::uint64_t>(edge_weight.weight));
  }
  return packed_values;
}

std::vector<CBC_SPDS::LocationEdgeWeight>
UnpackCBCSPDSLocationEdgeWeights(std::span<const std::uint64_t> packed_values)
{
  constexpr std::size_t packed_values_per_edge = 3;
  OpenSnLogicalErrorIf(packed_values.size() % packed_values_per_edge != 0,
                       "Malformed CBC location-edge-weight payload.");

  std::vector<CBC_SPDS::LocationEdgeWeight> edge_weights;
  edge_weights.reserve(packed_values.size() / packed_values_per_edge);
  for (std::size_t i = 0; i < packed_values.size(); i += packed_values_per_edge)
  {
    OpenSnLogicalErrorIf(
      packed_values[i] > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) or
        packed_values[i + 1] > static_cast<std::uint64_t>(std::numeric_limits<int>::max()),
      "Malformed CBC location-edge-weight rank.");
    edge_weights.push_back({static_cast<int>(packed_values[i]),
                            static_cast<int>(packed_values[i + 1]),
                            std::bit_cast<double>(packed_values[i + 2])});
  }
  return edge_weights;
}

void
AccumulateCBCGlobalEdgeWeights(const std::vector<std::shared_ptr<CBC_SPDS>>& spds_list)
{
  const int comm_size = opensn::mpi_comm.size();
  for (size_t spds_ordinal = 0; spds_ordinal < spds_list.size(); ++spds_ordinal)
  {
    const int owner = GetSweepGraphOwner(spds_ordinal, spds_list.size());
    auto local_values =
      PackCBCSPDSLocationEdgeWeights(spds_list[spds_ordinal]->ComputeLocalLocationEdgeWeights());
    const auto& dependencies = spds_list[spds_ordinal]->GetLocationDependencies();
    local_values.insert(local_values.end(), dependencies.begin(), dependencies.end());
    local_values.push_back(dependencies.size());
    OpenSnLogicalErrorIf(local_values.size() > static_cast<size_t>(std::numeric_limits<int>::max()),
                         "CBC location-edge-weight payload exceeds the MPI count limit.");

    const auto local_size = static_cast<int>(local_values.size());
    std::vector<int> receive_counts(comm_size, 0);
    mpi_comm.all_gather(local_size, receive_counts);

    std::vector<int> displacements(comm_size, 0);
    int total_size = 0;
    for (int rank = 0; rank < comm_size; ++rank)
    {
      displacements[rank] = total_size;
      OpenSnLogicalErrorIf(receive_counts[rank] > std::numeric_limits<int>::max() - total_size,
                           "CBC location-edge-weight gather exceeds the MPI count limit.");
      total_size += receive_counts[rank];
    }

    std::vector<std::uint64_t> global_values;
    mpi_comm.gather(local_values, global_values, receive_counts, displacements, owner);
    if (opensn::mpi_comm.rank() == owner)
    {
      std::vector<std::vector<int>> global_dependencies(comm_size);
      std::vector<CBC_SPDS::LocationEdgeWeight> global_weights;
      for (int rank = 0; rank < comm_size; ++rank)
      {
        const auto values = std::span<const std::uint64_t>(global_values)
                              .subspan(displacements[rank], receive_counts[rank]);
        OpenSnLogicalErrorIf(values.empty() or values.back() >= values.size(),
                             "Malformed CBC location-dependency payload.");
        const auto num_dependencies = static_cast<size_t>(values.back());
        const auto num_weight_values = values.size() - num_dependencies - 1;
        for (const auto dependency : values.subspan(num_weight_values, num_dependencies))
        {
          OpenSnLogicalErrorIf(std::cmp_greater_equal(dependency, comm_size),
                               "Malformed CBC location-dependency rank.");
          global_dependencies[rank].push_back(static_cast<int>(dependency));
        }
        const auto weights = UnpackCBCSPDSLocationEdgeWeights(values.first(num_weight_values));
        global_weights.insert(global_weights.end(), weights.begin(), weights.end());
      }
      spds_list[spds_ordinal]->SetGlobalDependencies(std::move(global_dependencies));
      spds_list[spds_ordinal]->SetGlobalEdgeWeights(global_weights);
    }
  }
}

void
BuildOwnedCBCSweepFAS(const std::vector<std::shared_ptr<CBC_SPDS>>& spds_list)
{
  log.Log0Verbose1() << program_timer.GetTimeString() << " Build global sweep FAS for CBC SPDS.";
  for (size_t spds_ordinal = 0; spds_ordinal < spds_list.size(); ++spds_ordinal)
    if (opensn::mpi_comm.rank() == GetSweepGraphOwner(spds_ordinal, spds_list.size()))
      spds_list[spds_ordinal]->BuildGlobalSweepFAS();
}

std::vector<int>
GatherCBCSweepFAS(const std::vector<std::shared_ptr<CBC_SPDS>>& spds_list)
{
  log.Log0Verbose1() << program_timer.GetTimeString() << " Gather FAS for CBC SPDS.";
  std::vector<int> local_edges_to_remove;
  for (size_t spds_ordinal = 0; spds_ordinal < spds_list.size(); ++spds_ordinal)
  {
    if (opensn::mpi_comm.rank() != GetSweepGraphOwner(spds_ordinal, spds_list.size()))
      continue;

    const auto& edges = spds_list[spds_ordinal]->GetGlobalSweepFAS();
    OpenSnLogicalErrorIf(spds_ordinal > static_cast<size_t>(std::numeric_limits<int>::max()) or
                           edges.size() > static_cast<size_t>(std::numeric_limits<int>::max()),
                         "CBC sweep FAS payload exceeds its integer serialization range.");
    local_edges_to_remove.push_back(static_cast<int>(spds_ordinal));
    local_edges_to_remove.push_back(static_cast<int>(edges.size()));
    local_edges_to_remove.insert(local_edges_to_remove.end(), edges.begin(), edges.end());
  }

  std::vector<int> receive_counts(opensn::mpi_comm.size(), 0);
  mpi_comm.all_gather(static_cast<int>(local_edges_to_remove.size()), receive_counts);

  std::vector<int> displacements(opensn::mpi_comm.size(), 0);
  int total_size = 0;
  for (size_t rank = 0; rank < receive_counts.size(); ++rank)
  {
    displacements[rank] = total_size;
    OpenSnLogicalErrorIf(receive_counts[rank] > std::numeric_limits<int>::max() - total_size,
                         "CBC sweep FAS gather exceeds the MPI count limit.");
    total_size += receive_counts[rank];
  }

  std::vector<int> global_edges_to_remove(total_size);
  mpi_comm.all_gather(local_edges_to_remove, global_edges_to_remove, receive_counts, displacements);
  return global_edges_to_remove;
}

void
ApplyCBCSweepFAS(const std::vector<std::shared_ptr<CBC_SPDS>>& spds_list,
                 const std::vector<int>& global_edges_to_remove)
{
  size_t offset = 0;
  while (offset < global_edges_to_remove.size())
  {
    const auto spds_ordinal = static_cast<size_t>(global_edges_to_remove[offset++]);
    const auto num_edges = global_edges_to_remove[offset++];

    std::vector<int> edges(global_edges_to_remove.begin() + static_cast<std::ptrdiff_t>(offset),
                           global_edges_to_remove.begin() +
                             static_cast<std::ptrdiff_t>(offset + num_edges));
    offset += static_cast<size_t>(num_edges);
    spds_list[spds_ordinal]->SetGlobalSweepFAS(std::move(edges));
  }

  log.Log0Verbose1() << program_timer.GetTimeString() << " Apply global sweep FAS for CBC SPDS.";
  for (const auto& spds : spds_list)
    spds->ApplyGlobalSweepFAS();
}

void
BuildCBCGlobalSweepGraph(SweepRuntime& runtime)
{
  auto spds_list = GetCBCSPDSList(runtime);
  std::erase_if(spds_list, [](const auto& spds) { return not spds->AllowCycles(); });
  if (spds_list.empty())
    return;

  AccumulateCBCGlobalEdgeWeights(spds_list);
  BuildOwnedCBCSweepFAS(spds_list);
  ApplyCBCSweepFAS(spds_list, GatherCBCSweepFAS(spds_list));
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
    std::shared_ptr<CBC_SPDS> spds;
  };
  std::vector<WorkItem> work;
  for (const auto& quadrature : runtime.quadrature_order)
  {
    const auto& spds_list = runtime.quadrature_spds_map.at(quadrature);
    for (const auto& spds : spds_list)
      work.push_back({quadrature, std::static_pointer_cast<CBC_SPDS>(spds)});
  }
  if (work.empty())
    return;

  std::ranges::sort(work, {}, [](const WorkItem& item) { return item.spds->GetId(); });

  const size_t nthreads = std::min<size_t>(std::max(1U, opensn_num_threads), work.size());
  log.Log() << program_timer.GetTimeString() << " FLUDS construction: " << work.size()
            << " angles (" << nthreads << " thread(s)).";
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
  {
    BuildCBCSPDS(
      runtime, groupsets, grid, face_neighbor_info, quadrature_allow_cycles_map, use_gpus);
    if (not use_gpus)
      BuildCBCGlobalSweepGraph(runtime);
  }
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
