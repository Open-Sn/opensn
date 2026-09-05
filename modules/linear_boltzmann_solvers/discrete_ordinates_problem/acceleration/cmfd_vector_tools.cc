// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_vector_tools.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace opensn
{
namespace
{

unsigned int
NumCoarseGroups(const unsigned int num_groups, const unsigned int group_aggregation_size)
{
  return (num_groups + group_aggregation_size - 1) / group_aggregation_size;
}

unsigned int
CoarseGroupBegin(const unsigned int coarse_group, const unsigned int group_aggregation_size)
{
  return coarse_group * group_aggregation_size;
}

unsigned int
CoarseGroupEnd(const unsigned int coarse_group,
               const unsigned int num_groups,
               const unsigned int group_aggregation_size)
{
  return std::min(num_groups, (coarse_group + 1) * group_aggregation_size);
}

} // namespace

std::vector<double>
CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                       const unsigned int first_group,
                       const unsigned int num_groups,
                       const CMFDCoarseMesh& coarse_mesh,
                       const std::vector<double>& phi)
{
  return CMFDRestrictScalarFlux(do_problem, first_group, num_groups, 1, coarse_mesh, phi);
}

std::vector<double>
CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                       const unsigned int first_group,
                       const unsigned int num_groups,
                       const unsigned int group_aggregation_size,
                       const CMFDCoarseMesh& coarse_mesh,
                       const std::vector<double>& phi)
{
  const auto& transport_views = do_problem.GetCellTransportViews();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();

  OpenSnInvalidArgumentIf(phi.size() != do_problem.GetPhiNewLocal().size(),
                          "Input scalar flux vector size does not match the transport problem.");

  OpenSnInvalidArgumentIf(group_aggregation_size == 0,
                          "CMFD group aggregation size must be greater than zero.");

  const unsigned int num_coarse_groups = NumCoarseGroups(num_groups, group_aggregation_size);
  std::vector<double> coarse_phi(coarse_mesh.NumLocalCells() * num_coarse_groups, 0.0);

  std::map<int, std::vector<uint64_t>> pid_keys;
  std::map<int, std::vector<double>> pid_values;

  for (const auto& membership : coarse_mesh.LocalFineCellMemberships())
  {
    std::vector<double> group_integrals(num_coarse_groups, 0.0);
    const auto& fine_cell = do_problem.GetGrid()->GetGlobalCell(membership.fine_cell_id);
    OpenSnInvalidArgumentIf(fine_cell.partition_id != opensn::mpi_comm.rank(),
                            "CMFD restriction fine-cell membership is not locally owned.");

    const auto& transport_view = transport_views[fine_cell.local_id];
    const auto& cell_matrices = unit_cell_matrices[fine_cell.local_id];

    for (int i = 0; i < transport_view.GetNumNodes(); ++i)
    {
      const auto node_volume = cell_matrices.intV_shapeI(i);
      const auto phi_map = transport_view.MapDOF(i, 0, first_group);
      for (unsigned int cg = 0; cg < num_coarse_groups; ++cg)
        for (unsigned int gsg = CoarseGroupBegin(cg, group_aggregation_size);
             gsg < CoarseGroupEnd(cg, num_groups, group_aggregation_size);
             ++gsg)
          group_integrals[cg] += phi[phi_map + gsg] * node_volume;
    }

    auto& keys = pid_keys[membership.coarse_cell_partition_id];
    auto& values = pid_values[membership.coarse_cell_partition_id];
    keys.push_back(membership.coarse_cell_id);
    values.push_back(fine_cell.volume);
    values.insert(values.end(), group_integrals.begin(), group_integrals.end());
  }

  const auto received_keys = MapAllToAll(pid_keys);
  const auto received_values = MapAllToAll(pid_values);

  std::vector<double> volumes(coarse_mesh.NumLocalCells(), 0.0);
  for (const auto& [pid, keys] : received_keys)
  {
    const auto values_it = received_values.find(pid);
    OpenSnLogicalErrorIf(values_it == received_values.end(),
                         "Missing CMFD restriction contribution values.");
    const auto& values = values_it->second;
    const auto record_size = static_cast<std::size_t>(num_coarse_groups) + 1;
    OpenSnLogicalErrorIf(values.size() != keys.size() * record_size,
                         "Invalid CMFD restriction contribution size.");

    for (std::size_t i = 0; i < keys.size(); ++i)
    {
      const auto& coarse_cell = coarse_mesh.LocalCellFromGlobalID(keys[i]);
      volumes[coarse_cell.local_id] += values[i * record_size];
      const auto coarse_offset = coarse_cell.local_id * num_coarse_groups;
      for (unsigned int cg = 0; cg < num_coarse_groups; ++cg)
        coarse_phi[coarse_offset + cg] += values[i * record_size + 1 + cg];
    }
  }

  for (const auto& coarse_cell : coarse_mesh.LocalCells())
  {
    const auto volume = volumes[coarse_cell.local_id];
    OpenSnLogicalErrorIf(volume <= 0.0, "CMFD coarse cell has non-positive volume.");
    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups;
    for (unsigned int cg = 0; cg < num_coarse_groups; ++cg)
      coarse_phi[coarse_offset + cg] /= volume;
  }

  return coarse_phi;
}

void
CMFDProlongateScalarFluxRatio(const DiscreteOrdinatesProblem& do_problem,
                              const unsigned int first_group,
                              const unsigned int num_groups,
                              const unsigned int group_aggregation_size,
                              const CMFDCoarseMesh& coarse_mesh,
                              const std::vector<double>& coarse_phi_old,
                              const std::vector<double>& coarse_phi_new,
                              std::vector<double>& phi)
{
  const auto& transport_views = do_problem.GetCellTransportViews();

  OpenSnInvalidArgumentIf(group_aggregation_size == 0,
                          "CMFD group aggregation size must be greater than zero.");
  const unsigned int num_coarse_groups = NumCoarseGroups(num_groups, group_aggregation_size);
  OpenSnInvalidArgumentIf(coarse_phi_old.size() != coarse_mesh.NumLocalCells() * num_coarse_groups,
                          "Old coarse scalar flux vector size mismatch.");
  OpenSnInvalidArgumentIf(coarse_phi_new.size() != coarse_mesh.NumLocalCells() * num_coarse_groups,
                          "New coarse scalar flux vector size mismatch.");
  OpenSnInvalidArgumentIf(phi.size() != do_problem.GetPhiNewLocal().size(),
                          "Output scalar flux vector size does not match the transport problem.");

  std::map<int, std::set<uint64_t>> pid_request_sets;
  for (const auto& membership : coarse_mesh.LocalFineCellMemberships())
    pid_request_sets[membership.coarse_cell_partition_id].insert(membership.coarse_cell_id);

  std::map<int, std::vector<uint64_t>> pid_requests;
  for (const auto& [pid, request_set] : pid_request_sets)
    pid_requests[pid] = {request_set.begin(), request_set.end()};

  const auto received_requests = MapAllToAll(pid_requests);

  std::map<int, std::vector<uint64_t>> pid_response_keys;
  std::map<int, std::vector<double>> pid_response_values;
  for (const auto& [pid, coarse_cell_ids] : received_requests)
  {
    auto& keys = pid_response_keys[pid];
    auto& values = pid_response_values[pid];
    for (const auto coarse_cell_id : coarse_cell_ids)
    {
      const auto& coarse_cell = coarse_mesh.LocalCellFromGlobalID(coarse_cell_id);
      const auto coarse_offset = coarse_cell.local_id * num_coarse_groups;
      keys.push_back(coarse_cell_id);
      for (unsigned int cg = 0; cg < num_coarse_groups; ++cg)
      {
        const double denominator = coarse_phi_old[coarse_offset + cg];
        const double ratio =
          std::fabs(denominator) > 1.0e-14 ? coarse_phi_new[coarse_offset + cg] / denominator : 1.0;
        values.push_back(ratio);
      }
    }
  }

  const auto received_response_keys = MapAllToAll(pid_response_keys);
  const auto received_response_values = MapAllToAll(pid_response_values);

  std::map<uint64_t, std::vector<double>> coarse_cell_ratios;
  for (const auto& [pid, keys] : received_response_keys)
  {
    const auto values_it = received_response_values.find(pid);
    OpenSnLogicalErrorIf(values_it == received_response_values.end(),
                         "Missing CMFD prolongation ratio values.");
    const auto& values = values_it->second;
    OpenSnLogicalErrorIf(values.size() != keys.size() * num_coarse_groups,
                         "Invalid CMFD prolongation ratio size.");

    for (std::size_t i = 0; i < keys.size(); ++i)
    {
      auto& ratios = coarse_cell_ratios[keys[i]];
      ratios.resize(num_coarse_groups);
      for (unsigned int cg = 0; cg < num_coarse_groups; ++cg)
        ratios[cg] = values[i * num_coarse_groups + cg];
    }
  }

  for (const auto& membership : coarse_mesh.LocalFineCellMemberships())
  {
    const auto ratios_it = coarse_cell_ratios.find(membership.coarse_cell_id);
    OpenSnLogicalErrorIf(ratios_it == coarse_cell_ratios.end(),
                         "Missing CMFD prolongation ratios for local fine-cell membership.");
    const auto& ratios = ratios_it->second;

    const auto& fine_cell = do_problem.GetGrid()->GetGlobalCell(membership.fine_cell_id);
    OpenSnInvalidArgumentIf(fine_cell.partition_id != opensn::mpi_comm.rank(),
                            "CMFD prolongation fine-cell membership is not locally owned.");

    const auto& transport_view = transport_views[fine_cell.local_id];
    for (int i = 0; i < transport_view.GetNumNodes(); ++i)
    {
      const auto phi_map = transport_view.MapDOF(i, 0, first_group);
      for (unsigned int cg = 0; cg < num_coarse_groups; ++cg)
      {
        for (unsigned int gsg = CoarseGroupBegin(cg, group_aggregation_size);
             gsg < CoarseGroupEnd(cg, num_groups, group_aggregation_size);
             ++gsg)
          phi[phi_map + gsg] *= ratios[cg];
      }
    }
  }
}

} // namespace opensn
