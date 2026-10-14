// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"
#include "framework/logging/log.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/utils/caliper_scopes.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

namespace opensn
{
namespace
{

constexpr unsigned int PRECURSOR_NODE_LAYOUT_VERSION = 2;
constexpr unsigned int CELLWISE_RESTART_VERSION = 3;

enum class CellwiseField
{
  PHI,
  PRECURSOR
};

struct RestartMeshData
{
  std::vector<uint64_t> cell_ids;
  std::vector<uint64_t> num_cell_nodes;
  std::vector<double> nodes_x;
  std::vector<double> nodes_y;
  std::vector<double> nodes_z;
};

std::filesystem::path
RankRestartPath(const std::filesystem::path& current_rank_path, int target_rank)
{
  auto path = current_rank_path.string();
  const auto suffix = std::to_string(opensn::mpi_comm.rank()) + std::string(".restart.h5");
  OpenSnInvalidArgumentIf(not path.ends_with(suffix),
                          "Restart path `" + path +
                            "` does not follow the <stem><rank>.restart.h5 convention.");
  path.resize(path.size() - suffix.size());
  return path + std::to_string(target_rank) + ".restart.h5";
}

bool
ReadRestartMeshData(hid_t file_id, RestartMeshData& data)
{
  bool success = true;
  success &= H5ReadDataset1D<uint64_t>(file_id, "mesh/cell_ids", data.cell_ids);
  success &= H5ReadDataset1D<uint64_t>(file_id, "mesh/num_cell_nodes", data.num_cell_nodes);
  success &= H5ReadDataset1D<double>(file_id, "mesh/nodes_x", data.nodes_x);
  success &= H5ReadDataset1D<double>(file_id, "mesh/nodes_y", data.nodes_y);
  success &= H5ReadDataset1D<double>(file_id, "mesh/nodes_z", data.nodes_z);
  if (not success)
    return false;

  OpenSnInvalidArgumentIf(data.cell_ids.size() != data.num_cell_nodes.size(),
                          "Restart mesh cell ID and node-count datasets have different sizes.");
  const auto num_nodes =
    std::accumulate(data.num_cell_nodes.begin(), data.num_cell_nodes.end(), uint64_t{0});
  OpenSnInvalidArgumentIf(data.nodes_x.size() != num_nodes or data.nodes_y.size() != num_nodes or
                            data.nodes_z.size() != num_nodes,
                          "Restart mesh coordinate datasets do not match the cell-node counts.");
  return true;
}

bool
WriteRestartMeshData(hid_t file_id, const LBSProblem& problem)
{
  const auto& grid = problem.GetGrid();
  const auto& discretization = problem.GetSpatialDiscretization();
  std::vector<uint64_t> cell_ids;
  std::vector<uint64_t> num_cell_nodes;
  std::vector<double> nodes_x;
  std::vector<double> nodes_y;
  std::vector<double> nodes_z;
  cell_ids.reserve(grid->GetLocalCellCount());
  num_cell_nodes.reserve(grid->GetLocalCellCount());

  for (const auto& cell : grid->GetLocalCells())
  {
    cell_ids.push_back(cell->global_id);
    const auto num_nodes = discretization.GetCellNumNodes(*cell);
    num_cell_nodes.push_back(num_nodes);
    const auto& nodes = discretization.GetCellNodeLocations(*cell);
    OpenSnLogicalErrorIf(nodes.size() != num_nodes,
                         problem.GetName() + ": inconsistent local cell-node metadata.");
    for (const auto& node : nodes)
    {
      nodes_x.push_back(node.x);
      nodes_y.push_back(node.y);
      nodes_z.push_back(node.z);
    }
  }

  bool success = H5CreateGroup(file_id, "mesh");
  success &= H5WriteDataset1D<uint64_t>(file_id, "mesh/cell_ids", cell_ids);
  success &= H5WriteDataset1D<uint64_t>(file_id, "mesh/num_cell_nodes", num_cell_nodes);
  success &= H5WriteDataset1D<double>(file_id, "mesh/nodes_x", nodes_x);
  success &= H5WriteDataset1D<double>(file_id, "mesh/nodes_y", nodes_y);
  success &= H5WriteDataset1D<double>(file_id, "mesh/nodes_z", nodes_z);
  return success;
}

bool
WriteCellwiseDataset(hid_t file_id,
                     const std::string& dataset_name,
                     const LBSProblem& problem,
                     const std::vector<double>& source,
                     CellwiseField field)
{
  const auto& grid = problem.GetGrid();
  const auto& discretization = problem.GetSpatialDiscretization();
  const auto& uk_man = problem.GetUnknownManager();
  const auto stride = field == CellwiseField::PHI ? problem.GetNumMoments() * problem.GetNumGroups()
                                                  : problem.GetNumPrecursors();
  std::vector<double> values;
  values.reserve(source.size());

  for (const auto& cell : grid->GetLocalCells())
    for (size_t i = 0; i < discretization.GetCellNumNodes(*cell); ++i)
      for (size_t c = 0; c < stride; ++c)
      {
        size_t dof = 0;
        if (field == CellwiseField::PHI)
        {
          const auto m = static_cast<unsigned int>(c / problem.GetNumGroups());
          const auto g = static_cast<unsigned int>(c % problem.GetNumGroups());
          dof = discretization.MapDOFLocal(*cell, i, uk_man, m, g);
        }
        else
          dof = discretization.MapDOFLocal(*cell, i) * stride + c;
        OpenSnLogicalErrorIf(dof >= source.size(),
                             problem.GetName() + ": restart source DOF is out of range.");
        values.push_back(source[dof]);
      }

  OpenSnLogicalErrorIf(values.size() != source.size(),
                       problem.GetName() + ": restart source has an incompatible size.");
  return H5WriteDataset1D<double>(file_id, dataset_name, values);
}

bool
SameLocalPartition(const LBSProblem& problem, const RestartMeshData& data)
{
  const auto& grid = problem.GetGrid();
  const auto& discretization = problem.GetSpatialDiscretization();
  if (data.cell_ids.size() != grid->GetLocalCellCount())
    return false;

  size_t c = 0;
  size_t node_offset = 0;
  for (const auto& cell : grid->GetLocalCells())
  {
    if (data.cell_ids[c] != cell->global_id or
        data.num_cell_nodes[c] != discretization.GetCellNumNodes(*cell))
      return false;
    const auto& local_nodes = discretization.GetCellNodeLocations(*cell);
    for (size_t i = 0; i < local_nodes.size(); ++i)
    {
      const Vector3 file_node(data.nodes_x[node_offset + i],
                              data.nodes_y[node_offset + i],
                              data.nodes_z[node_offset + i]);
      const auto coordinate_scale =
        std::max({1.0, std::abs(file_node.x), std::abs(file_node.y), std::abs(file_node.z)});
      const auto tolerance = 64.0 * std::numeric_limits<double>::epsilon() * coordinate_scale;
      if ((local_nodes[i] - file_node).NormSquare() > tolerance * tolerance)
        return false;
    }
    node_offset += local_nodes.size();
    ++c;
  }
  return true;
}

bool
ReadCellwiseDataset(const std::filesystem::path& current_rank_path,
                    int writer_mpi_size,
                    bool same_partition,
                    const std::string& dataset_name,
                    LBSProblem& problem,
                    std::vector<double>& destination,
                    CellwiseField field)
{
  const auto& grid = problem.GetGrid();
  const auto& discretization = problem.GetSpatialDiscretization();
  const auto& uk_man = problem.GetUnknownManager();
  const auto stride = field == CellwiseField::PHI ? problem.GetNumMoments() * problem.GetNumGroups()
                                                  : problem.GetNumPrecursors();
  std::vector<bool> found(grid->GetLocalCellCount(), false);
  std::fill(destination.begin(), destination.end(), 0.0);

  const int first_writer_rank = same_partition ? opensn::mpi_comm.rank() : 0;
  const int end_writer_rank = same_partition ? first_writer_rank + 1 : writer_mpi_size;
  for (int writer_rank = first_writer_rank; writer_rank < end_writer_rank; ++writer_rank)
  {
    const auto path = RankRestartPath(current_rank_path, writer_rank);
    const H5FileHandle file(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    OpenSnInvalidArgumentIf(
      file.Id() < 0, problem.GetName() + ": failed to open restart file `" + path.string() + "`.");

    RestartMeshData data;
    bool success = ReadRestartMeshData(file.Id(), data);
    std::vector<double> values;
    success &= H5ReadDataset1D<double>(file.Id(), dataset_name, values);
    if (not success)
      return false;

    const auto num_nodes =
      std::accumulate(data.num_cell_nodes.begin(), data.num_cell_nodes.end(), uint64_t{0});
    OpenSnInvalidArgumentIf(values.size() != num_nodes * stride,
                            problem.GetName() + ": restart dataset `" + dataset_name +
                              "` has size " + std::to_string(values.size()) +
                              " but its cell metadata requires " +
                              std::to_string(num_nodes * stride) + ".");

    size_t node_offset = 0;
    size_t value_offset = 0;
    for (size_t c = 0; c < data.cell_ids.size(); ++c)
    {
      const auto cell_id = data.cell_ids[c];
      const auto num_nodes_in_cell = data.num_cell_nodes[c];
      if (grid->IsCellLocal(cell_id))
      {
        const auto& cell = grid->GetGlobalCell(cell_id);
        OpenSnInvalidArgumentIf(found.at(cell.local_id),
                                problem.GetName() + ": restart contains duplicate cell " +
                                  std::to_string(cell_id) + ".");
        const auto& local_nodes = discretization.GetCellNodeLocations(cell);
        OpenSnInvalidArgumentIf(local_nodes.size() != num_nodes_in_cell,
                                problem.GetName() + ": restart cell " + std::to_string(cell_id) +
                                  " has an incompatible node count.");

        std::vector<size_t> file_to_local_node(num_nodes_in_cell, 0);
        std::vector<bool> local_node_used(num_nodes_in_cell, false);
        for (size_t i = 0; i < num_nodes_in_cell; ++i)
        {
          const Vector3 file_node(data.nodes_x[node_offset + i],
                                  data.nodes_y[node_offset + i],
                                  data.nodes_z[node_offset + i]);
          size_t best_node = num_nodes_in_cell;
          double best_distance = std::numeric_limits<double>::max();
          for (size_t j = 0; j < local_nodes.size(); ++j)
          {
            if (local_node_used[j])
              continue;
            const auto distance = (local_nodes[j] - file_node).NormSquare();
            if (distance < best_distance)
            {
              best_distance = distance;
              best_node = j;
            }
          }
          const auto coordinate_scale =
            std::max({1.0, std::abs(file_node.x), std::abs(file_node.y), std::abs(file_node.z)});
          const auto tolerance = 64.0 * std::numeric_limits<double>::epsilon() * coordinate_scale;
          OpenSnInvalidArgumentIf(
            best_node == num_nodes_in_cell or best_distance > tolerance * tolerance,
            problem.GetName() + ": restart node coordinates do not match cell " +
              std::to_string(cell_id) + ".");
          file_to_local_node[i] = best_node;
          local_node_used[best_node] = true;
        }

        for (size_t i = 0; i < num_nodes_in_cell; ++i)
          for (size_t component = 0; component < stride; ++component)
          {
            size_t dof = 0;
            if (field == CellwiseField::PHI)
            {
              const auto m = static_cast<unsigned int>(component / problem.GetNumGroups());
              const auto g = static_cast<unsigned int>(component % problem.GetNumGroups());
              dof = discretization.MapDOFLocal(cell, file_to_local_node[i], uk_man, m, g);
            }
            else
              dof = discretization.MapDOFLocal(cell, file_to_local_node[i]) * stride + component;
            destination.at(dof) = values.at(value_offset + i * stride + component);
          }
        found[cell.local_id] = true;
      }
      node_offset += num_nodes_in_cell;
      value_offset += num_nodes_in_cell * stride;
    }
  }

  const auto missing = std::find(found.begin(), found.end(), false);
  if (missing != found.end())
  {
    const auto local_id = std::distance(found.begin(), missing);
    const auto& cell = *std::next(grid->GetLocalCells().begin(), local_id);
    throw std::invalid_argument(problem.GetName() + ": restart data is missing global cell " +
                                std::to_string(cell->global_id) + ".");
  }
  return true;
}

bool
ReadSizedDoubleVector(hid_t file_id,
                      const std::string& dataset_name,
                      std::vector<double>& destination,
                      size_t expected_size,
                      const std::string& problem_name)
{
  std::vector<double> values;
  const bool success = H5ReadDataset1D<double>(file_id, dataset_name, values);
  if (success)
  {
    OpenSnInvalidArgumentIf(values.size() != expected_size,
                            problem_name + ": restart dataset `" + dataset_name + "` has size " +
                              std::to_string(values.size()) + " but expected " +
                              std::to_string(expected_size) + ".");
    destination = std::move(values);
  }
  return success;
}

bool
ReadPrecursorVector(hid_t file_id,
                    std::vector<double>& destination,
                    size_t expected_size,
                    const LBSProblem& problem,
                    unsigned int restart_format_version,
                    bool allow_size_remap)
{
  std::vector<double> values;
  const bool success = H5ReadDataset1D<double>(file_id, "precursors_new", values);
  if (not success)
    return false;

  const bool node_layout = restart_format_version >= PRECURSOR_NODE_LAYOUT_VERSION;
  if (node_layout and values.size() == expected_size)
  {
    destination = std::move(values);
    return true;
  }

  // Precursors are stored per spatial node (index node * J + j). Files written before this layout
  // stored one value per cell and family; expand those to the nodes of each cell.
  const auto& grid = problem.GetGrid();
  const auto& discretization = problem.GetSpatialDiscretization();
  const size_t num_local_nodes = discretization.GetNumLocalNodes();
  const size_t num_local_cells = grid->GetLocalCellCount();
  const size_t num_layout_entities = node_layout ? num_local_nodes : num_local_cells;
  OpenSnInvalidArgumentIf(num_layout_entities == 0 or num_local_nodes == 0,
                          problem.GetName() +
                            ": cannot remap restart precursor data without local cells and nodes.");
  OpenSnInvalidArgumentIf(
    values.size() % num_layout_entities != 0 or expected_size % num_local_nodes != 0,
    problem.GetName() + ": restart dataset `precursors_new` cannot be remapped from size " +
      std::to_string(values.size()) + " to " + std::to_string(expected_size) + ".");

  const size_t old_stride = values.size() / num_layout_entities;
  const size_t new_stride = expected_size / num_local_nodes;
  OpenSnInvalidArgumentIf(not allow_size_remap and old_stride != new_stride,
                          problem.GetName() + ": restart dataset `precursors_new` has size " +
                            std::to_string(values.size()) + " but expected " +
                            std::to_string(expected_size) + ".");
  if (not node_layout)
    log.Log0Warning() << problem.GetName()
                      << ": restart precursors use the former cell-averaged layout; expanding "
                         "each cell value to the nodes of that cell.";
  const size_t copy_stride = std::min(old_stride, new_stride);

  std::vector<double> remapped(expected_size, 0.0);
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = discretization.GetCellMapping(*cell);
    for (size_t i = 0; i < cell_mapping.GetNumNodes(); ++i)
    {
      const auto node_id = discretization.MapDOFLocal(*cell, i);
      const size_t old_base = (node_layout ? node_id : cell->local_id) * old_stride;
      const size_t new_base = node_id * new_stride;
      for (size_t j = 0; j < copy_stride; ++j)
        remapped[new_base + j] = values[old_base + j];
    }
  }

  destination = std::move(remapped);
  return true;
}

} // namespace

bool
LBSProblem::ReadRestartData(const RestartDataHook& extra_reader,
                            const std::filesystem::path& read_path,
                            bool allow_transient_initialization_from_steady)
{
  CaliperRegionScope cali_restart_io_scope("RestartIO", CaliperRestartIOScopeDepth());

  const auto& fname = read_path.empty() ? GetOptions().restart.read_path : read_path;
  OpenSnInvalidArgumentIf(fname.empty(), GetName() + ": restart read path is empty.");

  const auto header_path = RankRestartPath(fname, 0);
  const H5FileHandle header_file(H5Fopen(header_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  bool success = (header_file.Id() >= 0);
  OpenSnInvalidArgumentIf(
    not success, GetName() + ": failed to open restart file `" + header_path.string() + "`.");

  unsigned int restart_format_version = 1;
  success &= H5ReadOptionalAttribute<unsigned int>(
    header_file.Id(), "restart_format_version", restart_format_version);
  OpenSnInvalidArgumentIf(
    restart_format_version == 0 or restart_format_version > CELLWISE_RESTART_VERSION,
    GetName() + ": unsupported restart format version " + std::to_string(restart_format_version) +
      ". The newest supported version is " + std::to_string(CELLWISE_RESTART_VERSION) + ".");

  if (restart_format_version < CELLWISE_RESTART_VERSION)
  {
    const H5FileHandle file(H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    success &= (file.Id() >= 0);
    OpenSnInvalidArgumentIf(file.Id() < 0,
                            GetName() + ": failed to open restart file `" + fname.string() + "`.");

    const size_t expected_phi_size = phi_old_local_.size();
    const size_t expected_precursor_size = precursor_new_local_.size();
    if (H5Aexists(file.Id(), "mpi_size") > 0)
    {
      int restart_mpi_size = 0;
      success &= H5ReadAttribute<int>(file.Id(), "mpi_size", restart_mpi_size);
      OpenSnInvalidArgumentIf(restart_mpi_size != opensn::mpi_comm.size(),
                              GetName() + ": legacy restart was written with " +
                                std::to_string(restart_mpi_size) + " MPI ranks but this run uses " +
                                std::to_string(opensn::mpi_comm.size()) + ".");
    }
    if (H5Aexists(file.Id(), "mpi_rank") > 0)
    {
      int restart_mpi_rank = -1;
      success &= H5ReadAttribute<int>(file.Id(), "mpi_rank", restart_mpi_rank);
      OpenSnInvalidArgumentIf(restart_mpi_rank != opensn::mpi_comm.rank(),
                              GetName() + ": legacy restart file rank metadata is " +
                                std::to_string(restart_mpi_rank) + " but this process rank is " +
                                std::to_string(opensn::mpi_comm.rank()) + ".");
    }

    success &=
      ReadSizedDoubleVector(file.Id(), "phi_old", phi_old_local_, expected_phi_size, GetName());
    if (H5Has(file.Id(), "phi_new"))
      success &=
        ReadSizedDoubleVector(file.Id(), "phi_new", phi_new_local_, expected_phi_size, GetName());
    else
      phi_new_local_ = phi_old_local_;

    if (H5Has(file.Id(), "precursors_new"))
      success &= ReadPrecursorVector(file.Id(),
                                     precursor_new_local_,
                                     expected_precursor_size,
                                     *this,
                                     restart_format_version,
                                     allow_transient_initialization_from_steady);

    double time = GetTime();
    double dt = GetTimeStep();
    double theta = GetTheta();
    bool adjoint = GetOptions().adjoint;
    success &= H5ReadOptionalAttribute<double>(file.Id(), "time", time);
    success &= H5ReadOptionalAttribute<double>(file.Id(), "dt", dt);
    success &= H5ReadOptionalAttribute<double>(file.Id(), "theta", theta);
    success &= H5ReadOptionalAttribute<bool>(file.Id(), "adjoint", adjoint);
    OpenSnInvalidArgumentIf(adjoint != GetOptions().adjoint,
                            GetName() + ": restart adjoint mode does not match the configured "
                                        "problem mode.");
    if (success)
    {
      SetTime(time);
      SetTimeStep(dt);
      SetTheta(theta);
    }
    success &= ReadProblemRestartData(file.Id(), allow_transient_initialization_from_steady, false);
    if (extra_reader)
      success &= extra_reader(file.Id());
  }
  else
  {
    int writer_mpi_size = 0;
    unsigned int file_num_moments = 0;
    unsigned int file_num_groups = 0;
    unsigned int file_num_precursors = 0;
    success &= H5ReadAttribute<int>(header_file.Id(), "mpi_size", writer_mpi_size);
    success &= H5ReadAttribute<unsigned int>(header_file.Id(), "num_moments", file_num_moments);
    success &= H5ReadAttribute<unsigned int>(header_file.Id(), "num_groups", file_num_groups);
    success &=
      H5ReadAttribute<unsigned int>(header_file.Id(), "num_precursors", file_num_precursors);
    OpenSnInvalidArgumentIf(writer_mpi_size <= 0,
                            GetName() + ": restart contains an invalid writer MPI size.");
    OpenSnInvalidArgumentIf(file_num_moments != GetNumMoments() or
                              file_num_groups != GetNumGroups(),
                            GetName() + ": restart moment/group layout is incompatible with the "
                                        "configured problem.");
    OpenSnInvalidArgumentIf(file_num_precursors != GetNumPrecursors() and
                              not allow_transient_initialization_from_steady,
                            GetName() + ": restart precursor layout is incompatible with the "
                                        "configured problem.");

    bool same_partition = writer_mpi_size == opensn::mpi_comm.size();
    std::filesystem::path local_writer_path = header_path;
    if (opensn::mpi_comm.rank() < writer_mpi_size)
    {
      local_writer_path = RankRestartPath(fname, opensn::mpi_comm.rank());
      const H5FileHandle local_file(
        H5Fopen(local_writer_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      OpenSnInvalidArgumentIf(local_file.Id() < 0,
                              GetName() + ": failed to open restart file `" +
                                local_writer_path.string() + "`.");
      RestartMeshData local_data;
      success &= ReadRestartMeshData(local_file.Id(), local_data);
      same_partition &= SameLocalPartition(*this, local_data);
    }
    else
      same_partition = false;
    opensn::mpi_comm.all_reduce(same_partition, mpi::op::logical_and<bool>());

    success &= ReadCellwiseDataset(
      fname, writer_mpi_size, same_partition, "phi_old", *this, phi_old_local_, CellwiseField::PHI);
    if (H5Has(header_file.Id(), "phi_new"))
      success &= ReadCellwiseDataset(fname,
                                     writer_mpi_size,
                                     same_partition,
                                     "phi_new",
                                     *this,
                                     phi_new_local_,
                                     CellwiseField::PHI);
    else
      phi_new_local_ = phi_old_local_;

    if (H5Has(header_file.Id(), "precursors_new"))
      success &= ReadCellwiseDataset(fname,
                                     writer_mpi_size,
                                     same_partition,
                                     "precursors_new",
                                     *this,
                                     precursor_new_local_,
                                     CellwiseField::PRECURSOR);

    double time = GetTime();
    double dt = GetTimeStep();
    double theta = GetTheta();
    bool adjoint = GetOptions().adjoint;
    success &= H5ReadOptionalAttribute<double>(header_file.Id(), "time", time);
    success &= H5ReadOptionalAttribute<double>(header_file.Id(), "dt", dt);
    success &= H5ReadOptionalAttribute<double>(header_file.Id(), "theta", theta);
    success &= H5ReadOptionalAttribute<bool>(header_file.Id(), "adjoint", adjoint);
    OpenSnInvalidArgumentIf(adjoint != GetOptions().adjoint,
                            GetName() + ": restart adjoint mode does not match the configured "
                                        "problem mode.");
    if (success)
    {
      SetTime(time);
      SetTimeStep(dt);
      SetTheta(theta);
    }

    const H5FileHandle local_writer_file(
      H5Fopen(local_writer_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    OpenSnInvalidArgumentIf(local_writer_file.Id() < 0,
                            GetName() + ": failed to open restart file `" +
                              local_writer_path.string() + "`.");
    success &= ReadProblemRestartData(
      local_writer_file.Id(), allow_transient_initialization_from_steady, not same_partition);
    if (extra_reader)
      success &= extra_reader(header_file.Id());
  }

  if (success)
    log.Log() << "Successfully read restart data." << std::endl;
  else
    log.Log0Warning() << "Failed to read restart data from \"" << fname
                      << "\". Proceeding with the solver's default initial state instead of "
                         "the requested restart."
                      << std::endl;

  return success;
}

bool
LBSProblem::WriteRestartData(const RestartDataHook& extra_writer)
{
  CaliperRegionScope cali_restart_io_scope("RestartIO", CaliperRestartIOScopeDepth());

  const auto& fname = GetOptions().restart.write_path;
  OpenSnInvalidArgumentIf(fname.empty(), GetName() + ": restart write path is empty.");

  const H5FileHandle file(H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  bool success = (file.Id() >= 0);
  if (file.Id() >= 0)
  {
    success &= H5CreateAttribute<unsigned int>(
      file.Id(), "restart_format_version", CELLWISE_RESTART_VERSION);
    success &= H5CreateAttribute<int>(file.Id(), "mpi_size", opensn::mpi_comm.size());
    success &= H5CreateAttribute<int>(file.Id(), "mpi_rank", opensn::mpi_comm.rank());
    success &= H5CreateAttribute<unsigned int>(file.Id(), "num_moments", GetNumMoments());
    success &= H5CreateAttribute<unsigned int>(file.Id(), "num_groups", GetNumGroups());
    success &= H5CreateAttribute<unsigned int>(file.Id(), "num_precursors", GetNumPrecursors());
    success &= H5CreateAttribute<double>(file.Id(), "time", GetTime());
    success &= H5CreateAttribute<double>(file.Id(), "dt", GetTimeStep());
    success &= H5CreateAttribute<double>(file.Id(), "theta", GetTheta());
    success &= H5CreateAttribute<bool>(file.Id(), "adjoint", GetOptions().adjoint);
    success &= WriteRestartMeshData(file.Id(), *this);
    success &=
      WriteCellwiseDataset(file.Id(), "phi_old", *this, GetPhiOldLocal(), CellwiseField::PHI);
    success &=
      WriteCellwiseDataset(file.Id(), "phi_new", *this, GetPhiNewLocal(), CellwiseField::PHI);

    const auto& precursors_new_local = GetPrecursorsNewLocal();
    if (not precursors_new_local.empty())
      success &= WriteCellwiseDataset(
        file.Id(), "precursors_new", *this, precursors_new_local, CellwiseField::PRECURSOR);

    success &= WriteProblemRestartData(file.Id());
    if (extra_writer)
      success &= extra_writer(file.Id());
  }

  if (success)
  {
    UpdateRestartWriteTime();
    log.Log() << "Successfully wrote restart data." << std::endl;
  }
  else
    log.Log() << "Failed to write restart data." << std::endl;

  return success;
}

bool
LBSSolverIO::ReadRestartData(LBSProblem& lbs_problem,
                             const std::function<bool(hid_t)>& extra_reader,
                             const std::filesystem::path& read_path,
                             bool allow_transient_initialization_from_steady)
{
  return lbs_problem.ReadRestartData(
    extra_reader, read_path, allow_transient_initialization_from_steady);
}

bool
LBSSolverIO::WriteRestartData(LBSProblem& lbs_problem,
                              const std::function<bool(hid_t)>& extra_writer)
{
  return lbs_problem.WriteRestartData(extra_writer);
}

} // namespace opensn
