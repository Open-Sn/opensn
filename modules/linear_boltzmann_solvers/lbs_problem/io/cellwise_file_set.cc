// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/io/cellwise_file_set.h"
#include "framework/data_types/vector3.h"
#include "framework/logging/log.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/mpi/record_exchange.h"
#include "framework/utils/collective_error.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <unordered_set>

namespace opensn
{

namespace
{

/// The "mesh" group of one cell-wise file.
struct FileMesh
{
  std::vector<std::uint64_t> cell_ids;
  /// Start of each cell's nodes and the end of the last: the number of cells plus one entries once
  /// read, and empty before, so that an empty FileMesh allocates nothing.
  std::vector<std::uint64_t> node_offset;
  std::vector<double> nodes_x, nodes_y, nodes_z;

  std::uint64_t GetNumCellNodes(const std::size_t cell) const
  {
    return node_offset[cell + 1] - node_offset[cell];
  }

  /// Release the storage without allocating, leaving no cells and an empty `node_offset`.
  void Release() noexcept
  {
    std::vector<std::uint64_t>().swap(cell_ids);
    std::vector<std::uint64_t>().swap(node_offset);
    std::vector<double>().swap(nodes_x);
    std::vector<double>().swap(nodes_y);
    std::vector<double>().swap(nodes_z);
  }
};

/**
 * Read and check the "mesh" group of `file`, which may hold at most the `max_cells` cells and
 * `max_nodes` nodes of the whole current mesh. Throws on this rank.
 */
FileMesh
ReadFileMesh(const hid_t file,
             const std::string& file_name,
             const std::uint64_t max_cells,
             const std::uint64_t max_nodes)
{
  std::uint64_t num_cells = 0;
  std::uint64_t num_nodes = 0;
  OpenSnLogicalErrorIf(not H5ReadAttribute(file, "mesh/num_local_cells", num_cells) or
                         not H5ReadAttribute(file, "mesh/num_local_nodes", num_nodes),
                       "Failed to read the mesh metadata from " + file_name + ".");
  // The counts come from the file, so they are bounded by the mesh before anything is allocated.
  OpenSnLogicalErrorIf(num_cells > max_cells or num_nodes > max_nodes,
                       file_name + " records " + std::to_string(num_cells) + " cells and " +
                         std::to_string(num_nodes) + " nodes, more than the mesh has (" +
                         std::to_string(max_cells) + " cells and " + std::to_string(max_nodes) +
                         " nodes).");

  // Each dataset's extent is checked against the metadata before its buffer is allocated.
  FileMesh mesh;
  std::vector<std::uint64_t> num_cell_nodes;
  OpenSnLogicalErrorIf(
    not H5ReadDataset1D<std::uint64_t>(file, "mesh/cell_ids", mesh.cell_ids, num_cells) or
      not H5ReadDataset1D<std::uint64_t>(file, "mesh/num_cell_nodes", num_cell_nodes, num_cells),
    "Failed to read cell data matching the metadata from " + file_name + ".");
  OpenSnLogicalErrorIf(
    not H5ReadDataset1D<double>(file, "mesh/nodes_x", mesh.nodes_x, num_nodes) or
      not H5ReadDataset1D<double>(file, "mesh/nodes_y", mesh.nodes_y, num_nodes) or
      not H5ReadDataset1D<double>(file, "mesh/nodes_z", mesh.nodes_z, num_nodes),
    "Failed to read node coordinates matching the metadata from " + file_name + ".");
  // A non-finite coordinate would make the matching tolerance infinite and match any node.
  const auto is_finite = [](const std::vector<double>& values)
  { return std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); }); };
  OpenSnLogicalErrorIf(not is_finite(mesh.nodes_x) or not is_finite(mesh.nodes_y) or
                         not is_finite(mesh.nodes_z),
                       "The node coordinates in " + file_name + " contain a non-finite value.");

  std::unordered_set<std::uint64_t> seen_cells;
  mesh.node_offset.reserve(num_cells + 1);
  mesh.node_offset.push_back(0);
  for (std::uint64_t c = 0; c < num_cells; ++c)
  {
    OpenSnLogicalErrorIf(not seen_cells.insert(mesh.cell_ids[c]).second,
                         "Cell " + std::to_string(mesh.cell_ids[c]) +
                           " appears more than once in " + file_name + ".");
    OpenSnLogicalErrorIf(num_cell_nodes[c] == 0 or
                           num_cell_nodes[c] > num_nodes - mesh.node_offset.back(),
                         "The cells in " + file_name + " have more nodes than its metadata.");
    mesh.node_offset.push_back(mesh.node_offset.back() + num_cell_nodes[c]);
  }
  OpenSnLogicalErrorIf(mesh.node_offset.back() != num_nodes,
                       "The cells in " + file_name + " have fewer nodes than its metadata.");
  return mesh;
}

/**
 * Match each of the `num_nodes` file nodes at `xyz` (three coordinates per node) to a distinct node
 * of `nodes` at its location. Entry i of the result is the node matching file node i. Returns
 * nullopt if the nodes cannot be matched unambiguously. The coordinates must be finite.
 *
 * A file written from the same mesh stores each node location exactly, so file nodes are first
 * matched to nodes at bitwise the same location. Nodes at the same location, as in a collapsed
 * cell, pair in their stored order, which is the cell's node order on that mesh. Each remaining
 * file node must then have exactly one remaining node within a tolerance of its location, and no
 * two may share one; otherwise the match is ambiguous and rejected. Nodes at distinct locations are
 * therefore matched the same way whatever their order.
 */
std::optional<std::vector<std::size_t>>
MatchCellNodes(const double* xyz, const std::size_t num_nodes, const std::vector<Vector3>& nodes)
{
  if (num_nodes != nodes.size())
    return std::nullopt;
  constexpr auto unmatched = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> node_map(num_nodes, unmatched);
  std::vector<unsigned char> taken(nodes.size(), 0);
  const auto file_node = [&](const std::size_t n)
  { return Vector3(xyz[3 * n], xyz[3 * n + 1], xyz[3 * n + 2]); };

  // Exact matches.
  for (std::size_t n = 0; n < num_nodes; ++n)
  {
    const auto location = file_node(n);
    for (std::size_t m = 0; m < nodes.size(); ++m)
      if (taken[m] == 0 and nodes[m].x == location.x and nodes[m].y == location.y and
          nodes[m].z == location.z)
      {
        node_map[n] = m;
        taken[m] = 1;
        break;
      }
  }

  // Each remaining file node must have exactly one remaining candidate within the tolerance. The
  // distance is not squared, so that coordinates far apart cannot overflow to infinity and compare
  // equal to the tolerance.
  std::vector<std::size_t> candidate(num_nodes, unmatched);
  for (std::size_t n = 0; n < num_nodes; ++n)
  {
    if (node_map[n] != unmatched)
      continue;
    const auto location = file_node(n);
    for (std::size_t m = 0; m < nodes.size(); ++m)
    {
      if (taken[m] != 0)
        continue;
      const auto coordinate_scale = std::max({1.0,
                                              std::abs(location.x),
                                              std::abs(location.y),
                                              std::abs(location.z),
                                              std::abs(nodes[m].x),
                                              std::abs(nodes[m].y),
                                              std::abs(nodes[m].z)});
      const auto tolerance = 64.0 * std::numeric_limits<double>::epsilon() * coordinate_scale;
      const auto difference = nodes[m] - location;
      if (std::hypot(difference.x, difference.y, difference.z) <= tolerance)
      {
        if (candidate[n] != unmatched)
          return std::nullopt;
        candidate[n] = m;
      }
    }
    if (candidate[n] == unmatched)
      return std::nullopt;
  }
  // No two remaining file nodes may share their candidate.
  for (std::size_t n = 0; n < num_nodes; ++n)
  {
    if (candidate[n] == unmatched)
      continue;
    if (taken[candidate[n]] != 0)
      return std::nullopt;
    taken[candidate[n]] = 1;
    node_map[n] = candidate[n];
  }
  return node_map;
}

/// Marks a local cell that has not arrived from any file.
const std::uint64_t no_file = std::numeric_limits<std::uint64_t>::max();

/**
 * The files of a cell-wise file set as streams for RunFileSetPasses: stream 0 sends each cell's
 * node locations, with header {global id, node count}, and stream d + 1 sends the values of
 * dataset d, with header {global id}. A receiver matches a cell's nodes once and scatters each
 * dataset through that match.
 */
class CellwiseFileSource
{
public:
  CellwiseFileSource(
    const std::string& file_base,
    const bool single_file,
    const std::uint64_t num_files,
    const bool unrecorded_count,
    const std::vector<CellwiseDataset>& datasets,
    const Mesh& grid,
    const SpatialDiscretization& discretization,
    const std::function<void(hid_t, const std::string&)>& validate_file,
    const std::function<void(std::size_t, const Cell&, std::size_t, const double*)>& scatter,
    const OwnerDirectory& directory,
    const int num_ranks,
    const std::uint64_t file_set_id,
    const std::uint64_t max_cells,
    const std::uint64_t max_nodes)
    : file_base_(file_base),
      single_file_(single_file),
      num_files_(num_files),
      unrecorded_count_(unrecorded_count),
      datasets_(datasets),
      grid_(grid),
      discretization_(discretization),
      validate_file_(validate_file),
      scatter_(scatter),
      directory_(directory),
      num_ranks_(static_cast<std::uint64_t>(num_ranks)),
      file_set_id_(file_set_id),
      max_cells_(max_cells),
      max_nodes_(max_nodes)
  {
  }

  /// Allocate the per-cell state. Separate from the constructor, which allocates nothing, so that
  /// a failure here can be raised on every rank. Rank-local.
  void Allocate()
  {
    const auto num_cells = grid_.GetLocalCellCount();
    cell_file_.assign(num_cells, no_file);
    node_maps_.resize(num_cells);
    received_.assign(datasets_.size(), std::vector<unsigned char>(num_cells, 0));
  }

  void OpenFile(const std::uint64_t file_id)
  {
    file_name_ = FileName(file_id);
    file_ = H5FileHandle(H5Fopen(file_name_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    if (file_.Id() < 0)
    {
      std::string message = "Failed to open " + file_name_ + ".";
      if (unrecorded_count_ and not single_file_)
        message += " " + FileName(0) +
                   " does not record the number of files in its set, so the set is taken to "
                   "have one file per rank of this run; a set written by fewer ranks must be read "
                   "on as many ranks.";
      throw std::logic_error(message);
    }

    if (H5Aexists(file_.Id(), "num_partitions") > 0)
    {
      std::uint64_t count = 0;
      OpenSnLogicalErrorIf(not H5ReadAttribute(file_.Id(), "num_partitions", count),
                           "Failed to read the number of files from " + file_name_ + ".");
      OpenSnLogicalErrorIf(unrecorded_count_,
                           file_name_ + " records the number of files in its set, but " +
                             FileName(0) + " does not.");
      OpenSnLogicalErrorIf(count != num_files_,
                           file_name_ + " records a set of " + std::to_string(count) +
                             " file(s), but the set has " + std::to_string(num_files_) + ".");
      // Every file of one write records the same identifier.
      std::uint64_t id = 0;
      OpenSnLogicalErrorIf(H5Aexists(file_.Id(), "file_set_id") <= 0 or
                             not H5ReadAttribute(file_.Id(), "file_set_id", id),
                           file_name_ + " records the number of files in its set but no file-set "
                                        "identifier.");
      OpenSnLogicalErrorIf(id != file_set_id_,
                           file_name_ + " does not belong to the same write as " + FileName(0) +
                             ".");
    }
    else
      OpenSnLogicalErrorIf(not unrecorded_count_,
                           file_name_ + " does not record the number of files in its set, but " +
                             FileName(0) + " does.");

    validate_file_(file_.Id(), file_name_);
    mesh_ = ReadFileMesh(file_.Id(), file_name_, max_cells_, max_nodes_);
  }

  std::string PrepareFile(const bool has_file)
  {
    // The keys are passed by reference, so a rank copies nothing before the collective lookup.
    const std::vector<std::uint64_t> no_cells;
    owners_ = directory_.FindOwners(has_file ? mesh_.cell_ids : no_cells);
    std::string error;
    auto local_error = CaptureError(
      [&]
      {
        for (std::size_t c = 0; c < owners_.size(); ++c)
          if (owners_[c] < 0)
          {
            error = "Cell " + std::to_string(mesh_.cell_ids[c]) + " in " + file_name_ +
                    " is not a cell of the current mesh.";
            return;
          }
      });
    if (not error.empty())
      return error;
    return local_error;
  }

  std::size_t GetHeaderWidth(const std::size_t stream) const { return stream == 0 ? 2 : 1; }

  void BeginStream(const std::size_t stream)
  {
    next_cell_ = 0;
    std::vector<double>().swap(values_);
    if (stream == 0)
      return;
    const auto& dataset = datasets_[stream - 1];
    const auto num_nodes = mesh_.node_offset.back();
    OpenSnLogicalErrorIf(dataset.stride != 0 and
                           num_nodes > std::numeric_limits<std::uint64_t>::max() / dataset.stride,
                         "The size of " + dataset.name + " in " + file_name_ + " overflows.");
    // The extent is checked against the nodes in the metadata before the values are allocated.
    OpenSnLogicalErrorIf(
      not H5ReadDataset1D<double>(file_.Id(), dataset.name, values_, num_nodes * dataset.stride),
      "Failed to read " + dataset.name + " from " + file_name_ +
        ", or its size does not match the nodes in its metadata.");
    OpenSnLogicalErrorIf(std::any_of(values_.begin(),
                                     values_.end(),
                                     [](double value) { return not std::isfinite(value); }),
                         dataset.name + " in " + file_name_ + " contains a non-finite value.");
  }

  bool IsStreamDone(std::size_t /*stream*/) const { return next_cell_ == mesh_.cell_ids.size(); }

  void FillChunk(const std::size_t stream, const std::uint64_t budget, RecordBatch& batch)
  {
    const auto stride = stream == 0 ? 3 : datasets_[stream - 1].stride;
    std::uint64_t bytes = batch.GetPackedBytes();
    for (bool first = true; next_cell_ < mesh_.cell_ids.size(); first = false)
    {
      const auto c = next_cell_;
      const auto num_nodes = mesh_.GetNumCellNodes(c);
      const auto record_bytes = batch.GetRecordBytes(num_nodes * stride);
      if (not first and bytes + record_bytes > budget)
        break;
      bytes += record_bytes;

      auto* header = batch.AddRecord(owners_[c]);
      header[0] = static_cast<std::int64_t>(mesh_.cell_ids[c]);
      if (stream == 0)
      {
        header[1] = static_cast<std::int64_t>(num_nodes);
        for (auto node = mesh_.node_offset[c]; node < mesh_.node_offset[c + 1]; ++node)
        {
          const std::array<double, 3> xyz = {
            mesh_.nodes_x[node], mesh_.nodes_y[node], mesh_.nodes_z[node]};
          batch.AppendPayload(xyz.data(), xyz.size());
        }
      }
      else
        batch.AppendPayload(values_.data() + mesh_.node_offset[c] * stride, num_nodes * stride);
      ++records_sent_;
      ++next_cell_;
    }
  }

  void Receive(const std::size_t stream, const std::uint64_t pass, const ReceivedRecords& records)
  {
    for (std::size_t r = 0; r < records.GetNumRecords(); ++r)
    {
      ++records_received_;
      const auto* header = records.GetHeader(r);
      const auto cell_id = static_cast<std::uint64_t>(header[0]);
      const auto file_id = pass * num_ranks_ + static_cast<std::uint64_t>(records.sources[r]);
      OpenSnLogicalErrorIf(not grid_.IsCellLocal(cell_id),
                           "Cell " + std::to_string(cell_id) + " of " + FileName(file_id) +
                             " was sent to a rank that does not own it.");
      const auto local_id = grid_.MapCellGlobalID2LocalID(cell_id);
      const auto& cell = grid_.GetGlobalCell(cell_id);

      if (stream == 0)
      {
        OpenSnLogicalErrorIf(cell_file_[local_id] != no_file,
                             "Cell " + std::to_string(cell_id) + " appears in both " +
                               FileName(cell_file_[local_id]) + " and " + FileName(file_id) + ".");
        const auto num_nodes = static_cast<std::size_t>(header[1]);
        OpenSnLogicalErrorIf(records.GetPayloadSize(r) != 3 * num_nodes,
                             "Received inconsistent node locations of cell " +
                               std::to_string(cell_id) + ".");
        auto node_map = MatchCellNodes(
          records.GetPayload(r), num_nodes, discretization_.GetCellNodeLocations(cell));
        OpenSnLogicalErrorIf(not node_map.has_value(),
                             "The nodes of cell " + std::to_string(cell_id) + " in " +
                               FileName(file_id) + " do not match the nodes of the mesh.");
        node_maps_[local_id] = std::move(*node_map);
        cell_file_[local_id] = file_id;
        continue;
      }

      const auto dataset = stream - 1;
      const auto stride = datasets_[dataset].stride;
      const auto& node_map = node_maps_[local_id];
      OpenSnLogicalErrorIf(cell_file_[local_id] != file_id or received_[dataset][local_id] != 0 or
                             records.GetPayloadSize(r) != node_map.size() * stride,
                           "Received inconsistent " + datasets_[dataset].name + " values of cell " +
                             std::to_string(cell_id) + ".");
      const double* values = records.GetPayload(r);
      for (std::size_t i = 0; i < node_map.size(); ++i)
        scatter_(dataset, cell, node_map[i], values + i * stride);
      received_[dataset][local_id] = 1;
    }
  }

  /// Release the open file and its data without allocating, since it runs from a destructor.
  void CloseFile() noexcept
  {
    file_ = H5FileHandle(H5I_INVALID_HID);
    mesh_.Release();
    std::vector<int>().swap(owners_);
    std::vector<double>().swap(values_);
  }

  /// A description of the first local cell missing from the file set, or an empty string.
  std::string FindMissingCell() const
  {
    for (std::size_t local_id = 0; local_id < cell_file_.size(); ++local_id)
    {
      bool complete = cell_file_[local_id] != no_file;
      for (const auto& received : received_)
        complete = complete and received[local_id] != 0;
      if (complete)
        continue;
      std::string message = "Cell " + std::to_string(grid_.GetLocalCells()[local_id]->global_id) +
                            " has no data in the " + std::to_string(num_files_) + " file(s) of " +
                            file_base_ + ".";
      if (unrecorded_count_ and single_file_)
        message += " A single file that does not record the number of files in its set must hold "
                   "every cell of the mesh.";
      else if (unrecorded_count_)
        message += " The files do not record the number of files in their set, so it is taken to "
                   "be the rank count; a set written by more ranks must be read on as many ranks.";
      return message;
    }
    return {};
  }

  std::uint64_t GetRecordsSent() const { return records_sent_; }
  std::uint64_t GetRecordsReceived() const { return records_received_; }

private:
  std::string FileName(const std::uint64_t file_id) const
  {
    return single_file_ ? file_base_ + ".h5" : PartitionFileName(file_base_, file_id);
  }

  const std::string& file_base_;
  const bool single_file_;
  const std::uint64_t num_files_;
  const bool unrecorded_count_;
  const std::vector<CellwiseDataset>& datasets_;
  const Mesh& grid_;
  const SpatialDiscretization& discretization_;
  const std::function<void(hid_t, const std::string&)>& validate_file_;
  const std::function<void(std::size_t, const Cell&, std::size_t, const double*)>& scatter_;
  const OwnerDirectory& directory_;
  const std::uint64_t num_ranks_;
  const std::uint64_t file_set_id_;
  const std::uint64_t max_cells_;
  const std::uint64_t max_nodes_;

  H5FileHandle file_{H5I_INVALID_HID};
  std::string file_name_;
  FileMesh mesh_;
  std::vector<int> owners_;
  std::vector<double> values_;
  std::size_t next_cell_ = 0;

  /// Per local cell: the file it arrived from, and the local node of each of its file nodes.
  std::vector<std::uint64_t> cell_file_;
  std::vector<std::vector<std::size_t>> node_maps_;
  /// Per dataset and local cell: 1 once the cell's values have arrived.
  std::vector<std::vector<unsigned char>> received_;
  std::uint64_t records_sent_ = 0;
  std::uint64_t records_received_ = 0;
};

} // namespace

std::uint64_t
MakeFileSetId(const mpi::Communicator& comm)
{
  std::uint64_t file_set_id = 0;
  ThrowIfAnyRankFailed(comm,
                       comm.rank() == 0 ? CaptureError(
                                            [&]
                                            {
                                              std::random_device device;
                                              file_set_id =
                                                (static_cast<std::uint64_t>(device()) << 32) ^
                                                static_cast<std::uint64_t>(device());
                                            })
                                        : std::string(),
                       "Making a file-set identifier");
  comm.broadcast(file_set_id, 0);
  return file_set_id;
}

void
WriteCellwiseFileMesh(const hid_t file,
                      const Mesh& grid,
                      const SpatialDiscretization& discretization,
                      const mpi::Communicator& comm,
                      const std::uint64_t file_set_id)
{
  const auto num_local_cells = grid.GetLocalCellCount();
  const auto num_local_nodes = discretization.GetNumLocalNodes();

  std::vector<std::uint64_t> cell_ids, num_cell_nodes;
  cell_ids.reserve(num_local_cells);
  num_cell_nodes.reserve(num_local_cells);
  std::vector<double> nodes_x, nodes_y, nodes_z;
  nodes_x.reserve(num_local_nodes);
  nodes_y.reserve(num_local_nodes);
  nodes_z.reserve(num_local_nodes);
  for (const auto& cell : grid.GetLocalCells())
  {
    cell_ids.push_back(cell->global_id);
    num_cell_nodes.push_back(discretization.GetCellNumNodes(*cell));
    for (const auto& node : discretization.GetCellNodeLocations(*cell))
    {
      nodes_x.push_back(node.x);
      nodes_y.push_back(node.y);
      nodes_z.push_back(node.z);
    }
  }

  // The number of files in the set, so that any number of ranks can read it. A file missing it or
  // its mesh would only be rejected when the set is read, so a failed write is an error now.
  const auto check = [](const bool written, const std::string& name)
  { OpenSnLogicalErrorIf(not written, "Failed to write " + name + " to a cell-wise file."); };
  check(H5CreateAttribute(file, "num_partitions", static_cast<std::uint64_t>(comm.size())),
        "num_partitions");
  // Identifies the files of this write, so that a set mixing files of different writes is rejected.
  check(H5CreateAttribute(file, "file_set_id", file_set_id), "file_set_id");
  check(H5CreateGroup(file, "mesh"), "the mesh group");
  check(H5CreateAttribute(file, "mesh/num_local_cells", num_local_cells), "mesh/num_local_cells");
  check(H5CreateAttribute(file, "mesh/num_local_nodes", num_local_nodes), "mesh/num_local_nodes");
  check(H5WriteDataset1D(file, "mesh/cell_ids", cell_ids), "mesh/cell_ids");
  check(H5WriteDataset1D(file, "mesh/num_cell_nodes", num_cell_nodes), "mesh/num_cell_nodes");
  check(H5WriteDataset1D(file, "mesh/nodes_x", nodes_x), "mesh/nodes_x");
  check(H5WriteDataset1D(file, "mesh/nodes_y", nodes_y), "mesh/nodes_y");
  check(H5WriteDataset1D(file, "mesh/nodes_z", nodes_z), "mesh/nodes_z");
}

void
SetProvenanceCheck::Add(const SolutionProvenance& provenance, const std::string& file_name)
{
  OpenSnLogicalErrorIf(local_.has_value() and *local_ != provenance,
                       "The solution provenance in " + file_name +
                         " differs from that of the other files of the set.");
  local_ = provenance;
}

SolutionProvenance
SetProvenanceCheck::Agree(const mpi::Communicator& comm, const std::string& operation) const
{
  // One code per recorded (role, orientation) pair and one for a file without the attributes; a
  // rank without a file does not constrain the result.
  constexpr int num_orientations = static_cast<int>(AdjointFluxOrientation::PHYSICAL) + 1;
  constexpr int untagged = (static_cast<int>(FluxSolutionRole::ADJOINT) + 1) * num_orientations;
  constexpr int no_file_min = untagged + 1;
  constexpr int no_file_max = -1;
  std::optional<int> code;
  if (local_.has_value())
    code = local_->tagged ? static_cast<int>(local_->role) * num_orientations +
                              static_cast<int>(local_->orientation)
                          : untagged;
  int min_code = code.value_or(no_file_min);
  int max_code = code.value_or(no_file_max);
  comm.all_reduce(min_code, mpi::op::min<int>());
  comm.all_reduce(max_code, mpi::op::max<int>());
  if (min_code != max_code)
    throw std::runtime_error(operation +
                             " failed: the files of the set record different solution provenance: "
                             "solution roles, adjoint orientations, or whether they record them.");
  if (max_code == no_file_max or max_code == untagged)
    return {};
  return {true,
          static_cast<FluxSolutionRole>(max_code / num_orientations),
          static_cast<AdjointFluxOrientation>(max_code % num_orientations)};
}

SolutionProvenance
CheckSolutionProvenance(const hid_t file,
                        const std::string& file_name,
                        const FluxSolutionRole expected_role,
                        const AdjointFluxOrientation expected_orientation,
                        const std::optional<bool> problem_is_adjoint)
{
  auto file_role = static_cast<unsigned int>(FluxSolutionRole::UNSPECIFIED);
  auto file_orientation = static_cast<unsigned int>(AdjointFluxOrientation::UNSPECIFIED);
  const bool has_role = H5Aexists(file, "solution_role") > 0;
  const bool has_orientation = H5Aexists(file, "adjoint_orientation") > 0;
  if (has_role)
    OpenSnLogicalErrorIf(not H5ReadAttribute(file, "solution_role", file_role),
                         "Failed to read the solution role from " + file_name + ".");
  if (has_orientation)
    OpenSnLogicalErrorIf(not H5ReadAttribute(file, "adjoint_orientation", file_orientation),
                         "Failed to read the adjoint orientation from " + file_name + ".");
  OpenSnLogicalErrorIf(has_role != has_orientation,
                       "Incomplete solution provenance in " + file_name + ".");
  OpenSnLogicalErrorIf(file_role > static_cast<unsigned int>(FluxSolutionRole::ADJOINT),
                       "Invalid solution role in " + file_name + ".");
  OpenSnLogicalErrorIf(file_orientation >
                         static_cast<unsigned int>(AdjointFluxOrientation::PHYSICAL),
                       "Invalid adjoint orientation in " + file_name + ".");
  OpenSnLogicalErrorIf(
    has_role and
      ((file_role == static_cast<unsigned int>(FluxSolutionRole::ADJOINT)) !=
       (file_orientation != static_cast<unsigned int>(AdjointFluxOrientation::UNSPECIFIED))),
    "Inconsistent solution role and adjoint orientation in " + file_name + ".");
  OpenSnLogicalErrorIf(expected_role != FluxSolutionRole::UNSPECIFIED and
                         (not has_role or file_role != static_cast<unsigned int>(expected_role)),
                       "The solution role in " + file_name + " does not match the requested role.");
  OpenSnLogicalErrorIf(expected_orientation != AdjointFluxOrientation::UNSPECIFIED and
                         file_role == static_cast<unsigned int>(FluxSolutionRole::ADJOINT) and
                         file_orientation != static_cast<unsigned int>(expected_orientation),
                       "The adjoint orientation in " + file_name +
                         " does not match the requested orientation.");
  if (problem_is_adjoint.has_value() and has_role)
  {
    const auto problem_role =
      *problem_is_adjoint ? FluxSolutionRole::ADJOINT : FluxSolutionRole::FORWARD;
    OpenSnLogicalErrorIf(file_role != static_cast<unsigned int>(problem_role),
                         "The solution role in " + file_name + " does not match the problem mode.");
    OpenSnLogicalErrorIf(
      *problem_is_adjoint and
        file_orientation != static_cast<unsigned int>(AdjointFluxOrientation::SOLVER),
      "A physical adjoint flux in " + file_name +
        " cannot be loaded as problem state; use the restart interface so all solution "
        "vectors are restored to solver orientation together.");
  }
  return {has_role,
          static_cast<FluxSolutionRole>(file_role),
          static_cast<AdjointFluxOrientation>(file_orientation)};
}

void
ReadCellwiseFileSet(
  const std::string& file_base,
  const bool single_file,
  const std::vector<CellwiseDataset>& datasets,
  const Mesh& grid,
  const SpatialDiscretization& discretization,
  const std::function<void(hid_t, const std::string&)>& validate_file,
  const std::function<void(std::size_t, const Cell&, std::size_t, const double*)>& scatter,
  const mpi::Communicator& comm,
  const std::string& operation,
  const std::uint64_t max_send_bytes)
{
  OpenSnInvalidArgumentIf(max_send_bytes == 0, "The send budget must be positive.");

  // Rank 0 reads the number of files from the first file. Numbered files written before it was
  // recorded form a set of one file per rank of this run; a single file is always one file.
  // set_info holds the number of files, 1 if it is not recorded, and the file-set identifier.
  std::array<std::uint64_t, 3> set_info = {0, 0, 0};
  const auto read_set_info = [&]
  {
    const auto first_file = single_file ? file_base + ".h5" : PartitionFileName(file_base, 0);
    const H5FileHandle file(H5Fopen(first_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    OpenSnLogicalErrorIf(file.Id() < 0, "Failed to open " + first_file + ".");
    if (H5Aexists(file.Id(), "num_partitions") > 0)
    {
      OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), "num_partitions", set_info[0]) or
                             set_info[0] == 0,
                           "Failed to read a positive number of files from " + first_file + ".");
      // One file per writing rank, so no set has more files than an MPI run has ranks.
      OpenSnLogicalErrorIf(set_info[0] >
                             static_cast<std::uint64_t>(std::numeric_limits<int>::max()),
                           first_file + " records a set of " + std::to_string(set_info[0]) +
                             " files, more than an MPI run can have ranks.");
      OpenSnLogicalErrorIf(single_file and set_info[0] != 1,
                           first_file + " is one of a set of " + std::to_string(set_info[0]) +
                             " files; read the set instead of a single file.");
      OpenSnLogicalErrorIf(H5Aexists(file.Id(), "file_set_id") <= 0 or
                             not H5ReadAttribute(file.Id(), "file_set_id", set_info[2]),
                           first_file + " records the number of files in its set but no file-set "
                                        "identifier.");
    }
    else
      set_info[1] = 1;
  };
  // A failure on rank 0 is raised on every rank of `comm`.
  ThrowIfAnyRankFailed(
    comm, comm.rank() == 0 ? CaptureError(read_set_info) : std::string(), operation);
  comm.broadcast(set_info.data(), static_cast<int>(set_info.size()), 0);
  const bool unrecorded_count = set_info[1] != 0;
  const auto num_files =
    single_file ? 1 : (unrecorded_count ? static_cast<std::uint64_t>(comm.size()) : set_info[0]);

  // The whole mesh bounds the cells and nodes that any file may hold.
  const std::array<std::uint64_t, 2> local_size = {
    static_cast<std::uint64_t>(grid.GetLocalCellCount()),
    static_cast<std::uint64_t>(discretization.GetNumLocalNodes())};
  std::array<std::uint64_t, 2> mesh_size = {0, 0};
  comm.all_reduce(local_size.data(), 2, mesh_size.data(), mpi::op::sum<std::uint64_t>());

  // Rank-local setup runs inside CaptureError, and every rank checks the result before the next
  // collective, so a failure on one rank, such as an allocation, is raised on every rank.
  std::vector<std::uint64_t> local_cell_ids;
  ThrowIfAnyRankFailed(comm,
                       CaptureError(
                         [&]
                         {
                           local_cell_ids.reserve(grid.GetLocalCellCount());
                           for (const auto& cell : grid.GetLocalCells())
                             local_cell_ids.push_back(cell->global_id);
                         }),
                       operation);
  // OwnerDirectory throws std::logic_error on every rank together; raise it as documented.
  const OwnerDirectory directory = [&]
  {
    try
    {
      return OwnerDirectory(comm, local_cell_ids);
    }
    catch (const std::logic_error& exception)
    {
      throw std::runtime_error(operation + " failed: " + exception.what());
    }
  }();
  std::vector<std::uint64_t>().swap(local_cell_ids);

  CellwiseFileSource source(file_base,
                            single_file,
                            num_files,
                            unrecorded_count,
                            datasets,
                            grid,
                            discretization,
                            validate_file,
                            scatter,
                            directory,
                            comm.size(),
                            set_info[2],
                            mesh_size[0],
                            mesh_size[1]);
  ThrowIfAnyRankFailed(comm, CaptureError([&] { source.Allocate(); }), operation);
  // RunFileSetPasses reduces the budget on large rank counts to keep every exchange within the MPI
  // int count range.
  RunFileSetPasses<std::runtime_error>(
    comm, num_files, datasets.size() + 1, max_send_bytes, source, operation);

  std::string missing_cell;
  auto local_error = CaptureError([&] { missing_cell = source.FindMissingCell(); });
  ThrowIfAnyRankFailed(comm, missing_cell.empty() ? local_error : missing_cell, operation);
  std::array<std::uint64_t, 2> local_counts = {source.GetRecordsSent(),
                                               source.GetRecordsReceived()};
  std::array<std::uint64_t, 2> global_counts = {0, 0};
  comm.all_reduce(local_counts.data(), 2, global_counts.data(), mpi::op::sum<std::uint64_t>());
  // The counts are global, so every rank reaches the same decision.
  if (global_counts[0] != global_counts[1])
    throw std::runtime_error(operation + " failed: " + std::to_string(global_counts[0]) +
                             " records were sent, but " + std::to_string(global_counts[1]) +
                             " arrived.");
  // Best effort: a failed log must not keep this rank from the caller's next collective.
  try
  {
    log.Log0Verbose1() << operation << ": read " << num_files << " file(s) onto " << comm.size()
                       << " rank(s).";
  }
  catch (...) // NOLINT(bugprone-empty-catch)
  {
    // The log is dropped without building an error message, which could itself fail.
  }
}

} // namespace opensn
