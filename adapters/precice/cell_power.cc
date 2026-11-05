// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "adapters/precice/cell_power.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/geometry.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "framework/utils/collective_error.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace opensn::coupling
{
namespace
{
bool
IsFinite(const Vector3& value)
{
  return std::isfinite(value.x) and std::isfinite(value.y) and std::isfinite(value.z);
}
} // namespace

CellPowerData::CellPowerData(const unsigned int dimension, std::vector<OwnedCellPower> cells)
  : dimension_(dimension), cells_(std::move(cells))
{
  if (dimension_ < 1 or dimension_ > 3)
    throw std::invalid_argument("CellPowerData dimension must be 1, 2, or 3");
  std::unordered_set<std::uint64_t> cell_ids;
  cell_ids.reserve(cells_.size());
  for (const auto& cell : cells_)
  {
    if (not cell_ids.insert(cell.global_id).second)
      throw std::invalid_argument("CellPowerData contains a duplicate owned-cell global ID");
    if (not IsFinite(cell.centroid))
      throw std::invalid_argument("CellPowerData contains a non-finite cell centroid");
    if (not std::isfinite(cell.volume) or cell.volume <= 0.0)
      throw std::invalid_argument("CellPowerData contains a nonpositive or non-finite cell volume");
    if (not std::isfinite(cell.power))
      throw std::invalid_argument("CellPowerData contains non-finite cell power");
  }
}

double
CellPowerData::GetLocalPower() const
{
  double power = 0.0;
  for (const auto& cell : cells_)
    power += cell.power;
  return power;
}

void
CellPowerData::ConservativelyAverageByGroup(const std::vector<std::int64_t>& local_group_ids)
{
  int invalid_local_group_count = local_group_ids.size() != cells_.size() ? 1 : 0;
  mpi_comm.all_reduce(invalid_local_group_count, mpi::op::max<int>());
  if (invalid_local_group_count != 0)
    throw std::invalid_argument("local group count must equal the owned-cell count");

  // Gather only the distinct IDs of each rank so storage scales with the group count.
  std::vector<std::int64_t> local_unique_ids = local_group_ids;
  std::sort(local_unique_ids.begin(), local_unique_ids.end());
  local_unique_ids.erase(std::unique(local_unique_ids.begin(), local_unique_ids.end()),
                         local_unique_ids.end());
  std::size_t gathered_group_count = 0;
  mpi_comm.all_reduce(local_unique_ids.size(), gathered_group_count, mpi::op::sum<std::size_t>());
  if (gathered_group_count > static_cast<std::size_t>(std::numeric_limits<int>::max()) / 2)
    throw std::length_error("cell-power averaging group count exceeds the MPI count limit");
  if (gathered_group_count == 0)
    return;

  std::vector<std::int64_t> global_group_ids;
  mpi_comm.all_gather(local_unique_ids, global_group_ids);
  std::sort(global_group_ids.begin(), global_group_ids.end());
  global_group_ids.erase(std::unique(global_group_ids.begin(), global_group_ids.end()),
                         global_group_ids.end());

  std::vector<std::size_t> cell_group_indices(cells_.size());
  std::vector<double> local_totals(2 * global_group_ids.size(), 0.0);
  for (std::size_t cell_index = 0; cell_index < cells_.size(); ++cell_index)
  {
    const auto group = std::lower_bound(
      global_group_ids.begin(), global_group_ids.end(), local_group_ids[cell_index]);
    const auto group_index = static_cast<std::size_t>(group - global_group_ids.begin());
    cell_group_indices[cell_index] = group_index;
    local_totals[2 * group_index] += cells_[cell_index].power;
    local_totals[2 * group_index + 1] += cells_[cell_index].volume;
  }

  std::vector<double> global_totals(local_totals.size(), 0.0);
  mpi_comm.all_reduce(local_totals.data(),
                      static_cast<int>(local_totals.size()),
                      global_totals.data(),
                      mpi::op::sum<double>());
  for (std::size_t group = 0; group < global_group_ids.size(); ++group)
  {
    if (not std::isfinite(global_totals[2 * group]) or global_totals[2 * group] < 0.0)
      throw std::runtime_error("cell-power averaging group " +
                               std::to_string(global_group_ids[group]) +
                               " has negative or non-finite integrated power");
    if (not std::isfinite(global_totals[2 * group + 1]) or global_totals[2 * group + 1] <= 0.0)
      throw std::runtime_error("cell-power averaging group " +
                               std::to_string(global_group_ids[group]) +
                               " has invalid global volume");
  }

  for (std::size_t cell_index = 0; cell_index < cells_.size(); ++cell_index)
  {
    const auto group_index = cell_group_indices[cell_index];
    const double average_power_density =
      global_totals[2 * group_index] / global_totals[2 * group_index + 1];
    cells_[cell_index].power = average_power_density * cells_[cell_index].volume;
  }
}

CellPowerData
MakeCellPowerData(const FieldFunctionGridBased& field_function,
                  const std::size_t component,
                  const double source_rate,
                  const double joules_per_energy_unit,
                  const std::optional<double> transverse_measure,
                  const std::vector<unsigned int>& block_ids)
{
  if (not std::isfinite(source_rate) or source_rate < 0.0)
    throw std::invalid_argument("source_rate must be finite and nonnegative");
  if (not std::isfinite(joules_per_energy_unit) or joules_per_energy_unit <= 0.0)
    throw std::invalid_argument("joules_per_energy_unit must be finite and positive");

  const auto& unknown = field_function.GetUnknown();
  if (component >= unknown.num_components)
    throw std::out_of_range("field-function component is out of range");

  const auto& discretization = field_function.GetSpatialDiscretization();
  const auto grid = discretization.GetGrid();
  if (grid->GetCoordinateSystem() != CoordinateSystemType::CARTESIAN)
    throw std::invalid_argument("preCICE cell power requires a Cartesian mesh");
  const unsigned int dimension = grid->GetDimension();
  if (dimension < 1 or dimension > 3)
    throw std::invalid_argument("preCICE cell power requires a 1D, 2D, or 3D mesh");
  if (dimension == 3 and transverse_measure.has_value())
    throw std::invalid_argument("transverse_measure applies only to 1D and 2D meshes");
  if (dimension < 3 and not transverse_measure.has_value())
    throw std::invalid_argument(
      "transverse_measure is required for a 1D or 2D mesh: give the cross-sectional area of a "
      "1D mesh or the out-of-plane depth of a 2D mesh");
  const double transverse = transverse_measure.value_or(1.0);
  if (not std::isfinite(transverse) or transverse <= 0.0)
    throw std::invalid_argument("transverse_measure must be finite and positive");

  std::vector<unsigned int> selected_blocks = block_ids;
  std::sort(selected_blocks.begin(), selected_blocks.end());
  selected_blocks.erase(std::unique(selected_blocks.begin(), selected_blocks.end()),
                        selected_blocks.end());
  if (selected_blocks.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::length_error("block_ids exceeds the MPI count limit");
  if (not selected_blocks.empty())
  {
    std::vector<std::size_t> local_block_counts(selected_blocks.size(), 0);
    for (const auto& cell : grid->GetLocalCells())
    {
      const auto block =
        std::lower_bound(selected_blocks.begin(), selected_blocks.end(), cell->block_id);
      if (block != selected_blocks.end() and *block == cell->block_id)
        ++local_block_counts[block - selected_blocks.begin()];
    }
    std::vector<std::size_t> global_block_counts(selected_blocks.size(), 0);
    mpi_comm.all_reduce(local_block_counts.data(),
                        static_cast<int>(selected_blocks.size()),
                        global_block_counts.data(),
                        mpi::op::sum<std::size_t>());
    for (std::size_t b = 0; b < selected_blocks.size(); ++b)
      if (global_block_counts[b] == 0)
        throw std::invalid_argument("block_ids contains block " +
                                    std::to_string(selected_blocks[b]) +
                                    ", which has no cells in the mesh");
  }

  const auto& unknown_manager = field_function.GetUnknownManager();
  const auto& field_values = field_function.GetLocalFieldVector();
  const double conversion = source_rate * joules_per_energy_unit * transverse;

  // Failures below depend on rank-local data, so they are made collective before returning.
  CellPowerData cell_data(dimension, {});
  RunWithCollectiveErrors(
    [&]
    {
      std::vector<OwnedCellPower> cells;
      cells.reserve(grid->GetLocalCellCount());
      std::vector<double> node_values;
      for (const auto& cell : grid->GetLocalCells())
      {
        if (not selected_blocks.empty() and
            not std::binary_search(selected_blocks.begin(), selected_blocks.end(), cell->block_id))
          continue;

        const auto& cell_mapping = discretization.GetCellMapping(*cell);
        const std::size_t num_nodes = cell_mapping.GetNumNodes();
        node_values.assign(num_nodes, 0.0);
        for (std::size_t node = 0; node < num_nodes; ++node)
        {
          const auto dof = discretization.MapDOFLocal(
            *cell, node, unknown_manager, 0, static_cast<unsigned int>(component));
          if (dof >= field_values.size())
            throw std::logic_error("field-function local degree of freedom is out of range");
          node_values[node] = field_values[dof];
        }

        const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();
        double cell_integral = 0.0;
        double cell_measure = 0.0;
        for (const auto qp : fe_vol_data.GetQuadraturePointIndices())
        {
          double value = 0.0;
          for (std::size_t node = 0; node < num_nodes; ++node)
            value += fe_vol_data.ShapeValue(node, qp) * node_values[node];
          cell_integral += value * fe_vol_data.JxW(qp);
          cell_measure += fe_vol_data.JxW(qp);
        }

        cells.push_back(OwnedCellPower{cell->global_id,
                                       cell->block_id,
                                       cell->centroid,
                                       cell_measure * transverse,
                                       cell_integral * conversion});
      }
      cell_data = CellPowerData(dimension, std::move(cells));
    });

  return cell_data;
}

double
GetValidatedGlobalPower(const CellPowerData& cell_data)
{
  const auto& cells = cell_data.GetCells();
  std::size_t global_cell_count = 0;
  mpi_comm.all_reduce(cells.size(), global_cell_count, mpi::op::sum<std::size_t>());
  if (global_cell_count == 0)
    throw std::invalid_argument("the OpenSn preCICE coupling mesh has no cells");

  const auto dimension = static_cast<int>(cell_data.GetDimension());
  int minimum_dimension = 0;
  int maximum_dimension = 0;
  mpi_comm.all_reduce(dimension, minimum_dimension, mpi::op::min<int>());
  mpi_comm.all_reduce(dimension, maximum_dimension, mpi::op::max<int>());
  if (minimum_dimension != maximum_dimension)
    throw std::invalid_argument("OpenSn ranks report different coupling-mesh dimensions");

  double local_minimum_power = std::numeric_limits<double>::infinity();
  for (const auto& cell : cells)
    local_minimum_power = std::min(local_minimum_power, cell.power);
  double global_minimum_power = 0.0;
  mpi_comm.all_reduce(local_minimum_power, global_minimum_power, mpi::op::min<double>());
  if (global_minimum_power < 0.0)
    throw std::invalid_argument(
      "OpenSn cell-integrated deposition power must be nonnegative before transfer; "
      "use conservative averaging to remove finite-element undershoots");

  double global_power = 0.0;
  mpi_comm.all_reduce(cell_data.GetLocalPower(), global_power, mpi::op::sum<double>());
  if (not std::isfinite(global_power))
    throw std::invalid_argument("global OpenSn cell-integrated power must be finite");
  return global_power;
}

} // namespace opensn::coupling
