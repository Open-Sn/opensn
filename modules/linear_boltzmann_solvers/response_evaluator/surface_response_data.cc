// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/response_evaluator/surface_response_data.h"
#include "framework/logging/log.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace opensn
{

namespace
{

using FaceCentroidKey = DiscreteOrdinatesProblemIO::FaceCentroidKey;
using SurfaceAngularFlux = DiscreteOrdinatesProblemIO::SurfaceAngularFlux;

/// Tolerance for comparing quadrature values and direction cosines written by one problem.
const double same_problem_tolerance = 1.0e-12;

/**
 * Largest forward-to-adjoint direction distance, relative to the distance from the adjoint
 * direction to its nearest neighbor, at which the two are taken as the same direction. This is the
 * angular analog of matching face nodes within 1e-6 times the face size.
 */
const double direction_match_tolerance = 1.0e-6;

/// Throw on every rank if any rank recorded an error. Collective.
void
ThrowIfAnyRankFailed(const mpi::Communicator& comm,
                     const std::string& local_error,
                     const std::string& operation)
{
  int failed = local_error.empty() ? 0 : 1;
  comm.all_reduce(failed, mpi::op::max<int>());
  if (failed == 0)
    return;
  if (local_error.empty())
    throw std::logic_error(operation + " failed on another rank.");
  throw std::logic_error(operation + " failed on rank " + std::to_string(comm.rank()) + ": " +
                         local_error);
}

template <typename T>
void
BroadcastVector(const mpi::Communicator& comm, std::vector<T>& values, const int root)
{
  std::uint64_t size = values.size();
  comm.broadcast(size, root);
  values.resize(size);
  if (size > 0)
    comm.broadcast(values.data(), static_cast<int>(size), root);
}

/// Broadcast `strings` from `root` to every rank of `comm`. Collective.
void
BroadcastStrings(const mpi::Communicator& comm, std::vector<std::string>& strings, const int root)
{
  std::vector<std::uint64_t> lengths;
  std::vector<char> chars;
  for (const auto& string : strings)
  {
    lengths.push_back(string.size());
    chars.insert(chars.end(), string.begin(), string.end());
  }
  BroadcastVector(comm, lengths, root);
  BroadcastVector(comm, chars, root);
  strings.clear();
  std::uint64_t offset = 0;
  for (const auto length : lengths)
  {
    strings.emplace_back(chars.begin() + static_cast<std::ptrdiff_t>(offset),
                         chars.begin() + static_cast<std::ptrdiff_t>(offset + length));
    offset += length;
  }
}

std::string
PartitionFileName(const std::string& file_base, const std::uint64_t id)
{
  return file_base + std::to_string(id) + ".h5";
}

bool
SameQuadrature(const std::vector<SurfaceQuadrature>& a, const std::vector<SurfaceQuadrature>& b)
{
  const auto Close = [](const std::vector<double>& x, const std::vector<double>& y)
  {
    if (x.size() != y.size())
      return false;
    for (size_t i = 0; i < x.size(); ++i)
      if (std::abs(x[i] - y[i]) > same_problem_tolerance)
        return false;
    return true;
  };
  if (a.size() != b.size())
    return false;
  for (size_t gs = 0; gs < a.size(); ++gs)
    if (a[gs].num_dirs != b[gs].num_dirs or a[gs].num_groups != b[gs].num_groups or
        a[gs].first_group != b[gs].first_group or not Close(a[gs].omegas, b[gs].omegas) or
        not Close(a[gs].weights, b[gs].weights))
      return false;
  return true;
}

/// Attributes that every file of a surface angular-flux file set must share.
struct FileSetHeader
{
  std::uint64_t num_partitions = 0;
  std::uint64_t file_set_id = 0;
  bool adjoint = false;
  std::array<double, 3> centroid_origin = {0.0, 0.0, 0.0};
  double centroid_spacing = 0.0;
  /// Quadrature, group count, and first group of each groupset.
  std::vector<SurfaceQuadrature> quadratures;
  /// Per surface, 1 if the tag is one side of an interior surface and 0 for a boundary.
  std::vector<std::uint64_t> interior;
};

/**
 * Return the header of one file read by `DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes`,
 * whose entries are ordered by groupset, then by surface.
 */
FileSetHeader
GetFileHeader(const std::vector<SurfaceAngularFlux>& file_fluxes, const size_t num_surfaces)
{
  FileSetHeader header;
  const auto& first = file_fluxes.front();
  header.num_partitions = first.num_partitions;
  header.file_set_id = first.file_set_id;
  header.adjoint = first.adjoint;
  header.centroid_origin = first.mapping->centroid_origin;
  header.centroid_spacing = first.mapping->centroid_spacing;
  for (size_t i = 0; i < file_fluxes.size(); i += num_surfaces)
    header.quadratures.push_back({file_fluxes[i].num_directions,
                                  file_fluxes[i].num_groups,
                                  file_fluxes[i].first_group,
                                  file_fluxes[i].data.omega,
                                  file_fluxes[i].data.weights});
  for (size_t s = 0; s < num_surfaces; ++s)
    header.interior.push_back(file_fluxes[s].interior ? 1U : 0U);
  return header;
}

/// Return a description of how `header` differs from `reference`, or an empty string.
std::string
DescribeHeaderMismatch(const FileSetHeader& reference, const FileSetHeader& header)
{
  if (header.num_partitions != reference.num_partitions)
    return "it records " + std::to_string(header.num_partitions) +
           " partitions, but partition 0 records " + std::to_string(reference.num_partitions);
  if (header.file_set_id != reference.file_set_id)
    return "it was written by another WriteSurfaceAngularFluxes call than partition 0";
  if (header.adjoint != reference.adjoint)
    return "its adjoint flag differs from partition 0";
  if (header.centroid_origin != reference.centroid_origin or
      header.centroid_spacing != reference.centroid_spacing)
    return "its centroid origin or spacing differs from partition 0";
  if (not SameQuadrature(header.quadratures, reference.quadratures))
    return "its groupsets, quadratures, or group ranges differ from partition 0";
  if (header.interior != reference.interior)
    return "its surface kinds (boundary or interior side) differ from partition 0";
  return {};
}

/// Broadcast `header` from rank 0 to every rank of `comm`. Collective.
void
BroadcastFileSetHeader(const mpi::Communicator& comm, FileSetHeader& header)
{
  std::vector<std::uint64_t> ints;
  std::vector<double> reals;
  if (comm.rank() == 0)
  {
    ints = {header.num_partitions, header.file_set_id, header.adjoint ? 1U : 0U};
    reals = {header.centroid_origin[0],
             header.centroid_origin[1],
             header.centroid_origin[2],
             header.centroid_spacing};
    for (const auto& quadrature : header.quadratures)
    {
      ints.insert(ints.end(), {quadrature.num_dirs, quadrature.num_groups, quadrature.first_group});
      reals.insert(reals.end(), quadrature.omegas.begin(), quadrature.omegas.end());
      reals.insert(reals.end(), quadrature.weights.begin(), quadrature.weights.end());
    }
  }
  BroadcastVector(comm, ints, 0);
  BroadcastVector(comm, reals, 0);
  BroadcastVector(comm, header.interior, 0);

  header.num_partitions = ints[0];
  header.file_set_id = ints[1];
  header.adjoint = ints[2] != 0;
  header.centroid_origin = {reals[0], reals[1], reals[2]};
  header.centroid_spacing = reals[3];
  header.quadratures.resize((ints.size() - 3) / 3);
  auto real = reals.begin() + 4;
  for (size_t gs = 0; gs < header.quadratures.size(); ++gs)
  {
    auto& quadrature = header.quadratures[gs];
    quadrature.num_dirs = ints[3 + 3 * gs];
    quadrature.num_groups = ints[4 + 3 * gs];
    quadrature.first_group = ints[5 + 3 * gs];
    const auto num_dirs = static_cast<std::ptrdiff_t>(quadrature.num_dirs);
    quadrature.omegas.assign(real, real + 3 * num_dirs);
    real += 3 * num_dirs;
    quadrature.weights.assign(real, real + num_dirs);
    real += num_dirs;
  }
}

/// Faces read from this rank's files, in read order, before they are sent to their owners.
struct FaceChunk
{
  std::vector<FaceCentroidKey> keys;
  std::vector<std::uint64_t> num_nodes;
  std::vector<double> coords;
  std::vector<double> mass;
  /// Per groupset: one row of direction cosines per face.
  std::vector<std::vector<double>> mu;
  /// Per groupset: angular flux by face, node, direction, and group.
  std::vector<std::vector<double>> psi;
  std::uint64_t packed_bytes = 0;
};

/// Number of int64 header values per face: three key components and the node count.
const size_t face_header_size = 4;

/// Number of payload doubles of a face with `n` nodes.
std::uint64_t
FacePayloadSize(const std::uint64_t n, const std::vector<SurfaceQuadrature>& quadratures)
{
  std::uint64_t size = 3 * n + n * n;
  for (const auto& quadrature : quadratures)
    size += quadrature.num_dirs + n * quadrature.num_dirs * quadrature.num_groups;
  return size;
}

/**
 * The faces of one surface tag of one file read by
 * `DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes`, taken in bounded chunks.
 */
class TagFaces
{
public:
  TagFaces(const std::vector<SurfaceAngularFlux>& file_fluxes,
           const size_t surface_index,
           const size_t num_surfaces,
           const std::vector<SurfaceQuadrature>& quadratures)
    : quadratures_(quadratures)
  {
    for (size_t gs = 0; gs < quadratures.size(); ++gs)
      groupset_entries_.push_back(&file_fluxes[gs * num_surfaces + surface_index]);
    interior_ = groupset_entries_.front()->interior;
    const auto& mapping = *groupset_entries_.front()->mapping;
    keys_.resize(mapping.cell_ids.size());
    for (const auto& [key, face] : mapping.cell_map)
      keys_[face] = key;
    node_offset_.reserve(keys_.size() + 1);
    mass_offset_.reserve(keys_.size() + 1);
    for (const auto n : mapping.num_face_nodes)
    {
      node_offset_.push_back(node_offset_.back() + n);
      mass_offset_.push_back(mass_offset_.back() + n * n);
    }
  }

  bool IsDone() const { return next_face_ == keys_.size(); }

  /// Fill `chunk` with the next faces that fit in `budget` packed bytes, at least one face.
  void FillNextChunk(const std::uint64_t budget, FaceChunk& chunk)
  {
    std::uint64_t end = next_face_;
    while (end < keys_.size())
    {
      const auto bytes = sizeof(std::int64_t) * face_header_size +
                         sizeof(double) * FacePayloadSize(GetNumFaceNodes(end), quadratures_);
      if (end > next_face_ and chunk.packed_bytes + bytes > budget)
        break;
      chunk.packed_bytes += bytes;
      ++end;
    }
    if (end > next_face_)
      CopyFacesToChunk(end, chunk);
    next_face_ = end;
  }

private:
  std::uint64_t GetNumFaceNodes(const std::uint64_t face) const
  {
    return node_offset_[face + 1] - node_offset_[face];
  }

  /// Copy faces [next_face_, end) into `chunk`.
  void CopyFacesToChunk(const std::uint64_t end, FaceChunk& chunk) const
  {
    const auto& mapping = *groupset_entries_.front()->mapping;
    for (std::uint64_t f = next_face_; f < end; ++f)
    {
      chunk.keys.push_back(keys_[f]);
      chunk.num_nodes.push_back(GetNumFaceNodes(f));
      for (auto node = node_offset_[f]; node < node_offset_[f + 1]; ++node)
        chunk.coords.insert(chunk.coords.end(),
                            {mapping.nodes_x[node], mapping.nodes_y[node], mapping.nodes_z[node]});
    }

    const auto& mass = mapping.mass_matrix;
    chunk.mass.insert(chunk.mass.end(),
                      mass.begin() + static_cast<std::ptrdiff_t>(mass_offset_[next_face_]),
                      mass.begin() + static_cast<std::ptrdiff_t>(mass_offset_[end]));

    chunk.mu.resize(quadratures_.size());
    chunk.psi.resize(quadratures_.size());
    for (size_t gs = 0; gs < quadratures_.size(); ++gs)
    {
      const auto num_dirs = quadratures_[gs].num_dirs;
      const auto stride = num_dirs * quadratures_[gs].num_groups;
      const auto& data = groupset_entries_[gs]->data;
      // A boundary face carries its upwind values: those of the boundary condition for the
      // directions that enter the domain.
      const auto& face_psi = interior_ ? data.psi : data.upwind_psi;

      // The cosine of each direction with the face normal.
      const auto& omegas = quadratures_[gs].omegas;
      for (std::uint64_t f = next_face_; f < end; ++f)
      {
        const double* normal = &mapping.normals[3 * f];
        for (std::uint64_t d = 0; d < num_dirs; ++d)
          chunk.mu[gs].push_back(omegas[3 * d] * normal[0] + omegas[3 * d + 1] * normal[1] +
                                 omegas[3 * d + 2] * normal[2]);
      }
      chunk.psi[gs].insert(
        chunk.psi[gs].end(),
        face_psi.begin() + static_cast<std::ptrdiff_t>(node_offset_[next_face_] * stride),
        face_psi.begin() + static_cast<std::ptrdiff_t>(node_offset_[end] * stride));
    }
  }

  const std::vector<SurfaceQuadrature>& quadratures_;
  /// Entries of the tag, one per groupset.
  std::vector<const SurfaceAngularFlux*> groupset_entries_;
  bool interior_ = false;
  std::vector<FaceCentroidKey> keys_;
  std::vector<std::uint64_t> node_offset_ = {0};
  std::vector<std::uint64_t> mass_offset_ = {0};
  std::uint64_t next_face_ = 0;
};

std::vector<int>
MakeOffsets(const std::vector<int>& counts)
{
  std::vector<int> offsets(counts.size(), 0);
  for (size_t p = 1; p < counts.size(); ++p)
    offsets[p] = offsets[p - 1] + counts[p - 1];
  return offsets;
}

/**
 * Send the faces of `chunk` to their owners and append the received faces to `shard`. Collective.
 * Returns a description of a local unpacking error, or an empty string.
 */
std::string
ExchangeChunk(const mpi::Communicator& comm,
              const std::vector<SurfaceQuadrature>& quadratures,
              FaceChunk& chunk,
              SurfaceShard& shard)
{
  const auto num_ranks = static_cast<size_t>(comm.size());
  const auto num_faces = chunk.keys.size();

  // Count header and payload values per destination rank.
  std::vector<size_t> owner(num_faces);
  std::vector<std::uint64_t> send_counts(2 * num_ranks, 0);
  for (size_t f = 0; f < num_faces; ++f)
  {
    owner[f] = static_cast<size_t>(SurfaceFaceOwner(chunk.keys[f], comm.size()));
    send_counts[2 * owner[f]] += face_header_size;
    send_counts[2 * owner[f] + 1] += FacePayloadSize(chunk.num_nodes[f], quadratures);
  }
  std::vector<std::uint64_t> recv_counts(2 * num_ranks, 0);
  comm.all_to_all(send_counts.data(), 2, recv_counts.data(), 2);

  // The send budget bounds what each rank sends, but a rank receives from every rank, so check the
  // totals collectively before using int counts.
  std::uint64_t recv_headers = 0;
  std::uint64_t recv_payload = 0;
  for (size_t p = 0; p < num_ranks; ++p)
  {
    recv_headers += recv_counts[2 * p];
    recv_payload += recv_counts[2 * p + 1];
  }
  std::uint64_t send_payload = 0;
  for (size_t p = 0; p < num_ranks; ++p)
    send_payload += send_counts[2 * p + 1];
  const auto int_max = static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  ThrowIfAnyRankFailed(comm,
                       std::max({recv_headers, recv_payload, send_payload}) > int_max
                         ? "one exchange exceeds the MPI int count range of a rank"
                         : "",
                       "Surface angular-flux redistribution");

  std::vector<int> header_counts(num_ranks), payload_counts(num_ranks);
  std::vector<int> recv_header_counts(num_ranks), recv_payload_counts(num_ranks);
  for (size_t p = 0; p < num_ranks; ++p)
  {
    header_counts[p] = static_cast<int>(send_counts[2 * p]);
    payload_counts[p] = static_cast<int>(send_counts[2 * p + 1]);
    recv_header_counts[p] = static_cast<int>(recv_counts[2 * p]);
    recv_payload_counts[p] = static_cast<int>(recv_counts[2 * p + 1]);
  }
  const auto header_offsets = MakeOffsets(header_counts);
  const auto payload_offsets = MakeOffsets(payload_counts);

  // Pack in destination order.
  std::vector<std::int64_t> headers(num_faces * face_header_size);
  std::vector<double> payload(send_payload);
  {
    auto header_cursor = header_offsets;
    auto payload_cursor = payload_offsets;
    std::uint64_t node_offset = 0;
    std::uint64_t mass_offset = 0;
    std::vector<std::uint64_t> psi_offsets(quadratures.size(), 0);
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto n = chunk.num_nodes[f];
      auto* h = &headers[header_cursor[owner[f]]];
      h[0] = std::get<0>(chunk.keys[f]);
      h[1] = std::get<1>(chunk.keys[f]);
      h[2] = std::get<2>(chunk.keys[f]);
      h[3] = static_cast<std::int64_t>(n);
      header_cursor[owner[f]] += static_cast<int>(face_header_size);

      auto* out = &payload[payload_cursor[owner[f]]];
      out = std::copy_n(&chunk.coords[3 * node_offset], 3 * n, out);
      out = std::copy_n(&chunk.mass[mass_offset], n * n, out);
      for (size_t gs = 0; gs < quadratures.size(); ++gs)
      {
        const auto num_dirs = quadratures[gs].num_dirs;
        const auto size = n * num_dirs * quadratures[gs].num_groups;
        out = std::copy_n(&chunk.mu[gs][f * num_dirs], num_dirs, out);
        out = std::copy_n(&chunk.psi[gs][psi_offsets[gs]], size, out);
        psi_offsets[gs] += size;
      }
      payload_cursor[owner[f]] += static_cast<int>(FacePayloadSize(n, quadratures));
      node_offset += n;
      mass_offset += n * n;
    }
  }
  chunk = FaceChunk();

  std::vector<std::int64_t> recv_header_values;
  comm.all_to_all(headers,
                  header_counts,
                  header_offsets,
                  recv_header_values,
                  recv_header_counts,
                  MakeOffsets(recv_header_counts));
  headers = {};
  std::vector<double> recv_payload_values;
  comm.all_to_all(payload,
                  payload_counts,
                  payload_offsets,
                  recv_payload_values,
                  recv_payload_counts,
                  MakeOffsets(recv_payload_counts));
  payload = {};

  // Unpack in arrival order.
  const double* in = recv_payload_values.data();
  const double* const in_end = in + recv_payload_values.size();
  for (size_t h = 0; h < recv_header_values.size(); h += face_header_size)
  {
    const auto n = static_cast<std::uint64_t>(recv_header_values[h + 3]);
    if (n == 0 or static_cast<std::uint64_t>(in_end - in) < FacePayloadSize(n, quadratures))
      return "received an inconsistent surface face";
    const FaceCentroidKey key = {
      recv_header_values[h], recv_header_values[h + 1], recv_header_values[h + 2]};
    shard.index.emplace_back(key, shard.GetNumFaces());
    shard.coords.insert(shard.coords.end(), in, in + 3 * n);
    in += 3 * n;
    shard.mass.insert(shard.mass.end(), in, in + n * n);
    in += n * n;
    for (size_t gs = 0; gs < quadratures.size(); ++gs)
    {
      const auto num_dirs = quadratures[gs].num_dirs;
      const auto size = n * num_dirs * quadratures[gs].num_groups;
      auto& values = shard.groupsets[gs];
      values.mu.insert(values.mu.end(), in, in + num_dirs);
      in += num_dirs;
      values.psi.insert(values.psi.end(), in, in + size);
      in += size;
    }
    shard.node_offset.push_back(shard.node_offset.back() + n);
    shard.mass_offset.push_back(shard.mass_offset.back() + n * n);
  }
  if (in != in_end)
    return "received surface payload that does not match its face headers";
  return {};
}

} // namespace

std::optional<std::uint64_t>
SurfaceShard::Find(const DiscreteOrdinatesProblemIO::FaceCentroidKey& key) const
{
  const auto it =
    std::lower_bound(index.begin(),
                     index.end(),
                     key,
                     [](const auto& entry, const auto& value) { return entry.first < value; });
  if (it == index.end() or it->first != key)
    return std::nullopt;
  return it->second;
}

int
SurfaceFaceOwner(const DiscreteOrdinatesProblemIO::FaceCentroidKey& key, const int num_ranks)
{
  // splitmix64 finalizer: deterministic across ranks, builds, and platforms.
  const auto Mix = [](std::uint64_t x)
  {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  };
  auto hash = Mix(static_cast<std::uint64_t>(std::get<0>(key)));
  hash = Mix(hash ^ static_cast<std::uint64_t>(std::get<1>(key)));
  hash = Mix(hash ^ static_cast<std::uint64_t>(std::get<2>(key)));
  return static_cast<int>(hash % static_cast<std::uint64_t>(num_ranks));
}

DistributedSurfaceData
DistributeSurfaceAngularFluxes(const std::string& file_base,
                               const std::vector<std::string>& surfaces,
                               const mpi::Communicator& comm,
                               const std::uint64_t max_send_bytes)
{
  OpenSnInvalidArgumentIf(surfaces.empty(), "At least one surface name is required.");
  OpenSnInvalidArgumentIf(std::set<std::string>(surfaces.begin(), surfaces.end()).size() !=
                            surfaces.size(),
                          "Surface names must be unique.");
  OpenSnInvalidArgumentIf(max_send_bytes == 0, "The send budget must be positive.");
  const std::string operation = "Reading surface angular fluxes from " + file_base;
  const auto rank = comm.rank();
  const auto num_ranks = comm.size();

  // Rank 0 reads file 0, which resolves the surface names to stored tags and gives the attributes
  // of the file set. Every other file is read with the same tags.
  std::vector<SurfaceAngularFlux> file_fluxes;
  std::vector<std::string> surface_tags;
  FileSetHeader header;
  std::vector<char> root_error;
  if (rank == 0)
  {
    try
    {
      file_fluxes = DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(file_base, 0, surfaces);
      for (const auto& surf_flux : file_fluxes)
        if (surf_flux.groupset_id == 0)
          surface_tags.push_back(surf_flux.surface_name);
      header = GetFileHeader(file_fluxes, surface_tags.size());
    }
    catch (const std::exception& error)
    {
      const std::string message = error.what();
      root_error.assign(message.begin(), message.end());
    }
  }
  BroadcastVector(comm, root_error, 0);
  if (not root_error.empty())
    throw std::logic_error(operation +
                           " failed: " + std::string(root_error.begin(), root_error.end()));

  BroadcastStrings(comm, surface_tags, 0);
  BroadcastFileSetHeader(comm, header);
  const auto num_surfaces = surface_tags.size();

  // Surface responses pair both sides of an interior surface, so a side is read only with its
  // other side. The tags and kinds are identical on all ranks.
  for (size_t s = 0; s < num_surfaces; ++s)
  {
    if (header.interior[s] == 0)
      continue;
    const auto& tag = surface_tags[s];
    const auto other = tag.substr(0, tag.size() - 2) + (tag.ends_with("_u") ? "_d" : "_u");
    OpenSnInvalidArgumentIf(
      std::find(surface_tags.begin(), surface_tags.end(), other) == surface_tags.end(),
      std::string("Side ")
        .append(tag)
        .append(" of an interior surface was requested without its other side ")
        .append(other)
        .append(". Request the interior surface by its name, ")
        .append(tag.substr(0, tag.size() - 2))
        .append("."));
  }
  const auto num_partitions = header.num_partitions;
  const auto num_groupsets = header.quadratures.size();

  DistributedSurfaceData data;
  data.adjoint = header.adjoint;
  data.num_partitions = num_partitions;
  data.centroid_origin = header.centroid_origin;
  data.centroid_spacing = header.centroid_spacing;
  // Every file of the set must have this quadrature; see DescribeHeaderMismatch.
  data.groupsets = header.quadratures;
  for (size_t s = 0; s < num_surfaces; ++s)
  {
    auto& shard = data.surfaces[surface_tags[s]];
    shard.interior = header.interior[s] != 0;
    shard.groupsets.resize(num_groupsets);
  }

  // File k is read by rank k mod P, as its (k / P)-th file, and its faces are sent to their
  // owners before the next file is read. Every rank runs the loop max_files_per_rank times, so
  // ranks without a file in a pass still join its collectives.
  std::string error;
  std::vector<std::uint64_t> faces_sent(num_surfaces, 0);
  const auto num_ranks_u64 = static_cast<std::uint64_t>(num_ranks);
  const auto max_files_per_rank = (num_partitions + num_ranks_u64 - 1) / num_ranks_u64;
  // This rank's i-th file.
  for (std::uint64_t i = 0; i < max_files_per_rank; ++i)
  {
    const auto file_id = i * num_ranks_u64 + static_cast<std::uint64_t>(rank);
    const auto file_name = PartitionFileName(file_base, file_id);
    try
    {
      // Rank 0 read file 0 above.
      if (file_id > 0)
      {
        file_fluxes.clear();
        if (file_id < num_partitions)
        {
          file_fluxes =
            DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(file_base, file_id, surface_tags);
          const auto mismatch =
            DescribeHeaderMismatch(header, GetFileHeader(file_fluxes, num_surfaces));
          OpenSnLogicalErrorIf(not mismatch.empty(),
                               std::string(file_name)
                                 .append(" does not belong to the file set: ")
                                 .append(mismatch)
                                 .append("."));
        }
      }
    }
    catch (const std::exception& exception)
    {
      error = exception.what();
      file_fluxes.clear();
    }
    ThrowIfAnyRankFailed(comm, error, operation);

    for (size_t s = 0; s < num_surfaces; ++s)
    {
      auto& shard = data.surfaces[surface_tags[s]];
      std::optional<TagFaces> faces;
      if (error.empty() and not file_fluxes.empty())
      {
        try
        {
          faces.emplace(file_fluxes, s, num_surfaces, data.groupsets);
        }
        catch (const std::exception& exception)
        {
          error = exception.what();
        }
      }
      while (true)
      {
        FaceChunk chunk;
        if (error.empty() and faces.has_value() and not faces->IsDone())
        {
          try
          {
            faces->FillNextChunk(max_send_bytes, chunk);
          }
          catch (const std::exception& exception)
          {
            error = exception.what();
            chunk = FaceChunk();
          }
        }
        faces_sent[s] += chunk.keys.size();
        // Every rank takes part in every exchange until no rank has faces left.
        std::array<int, 2> local_flags = {error.empty() ? 0 : 1, chunk.keys.empty() ? 0 : 1};
        std::array<int, 2> global_flags = {0, 0};
        comm.all_reduce(local_flags.data(), 2, global_flags.data(), mpi::op::max<int>());
        if (global_flags[0] != 0)
          ThrowIfAnyRankFailed(comm, error, operation);
        if (global_flags[1] == 0)
          break;
        error = ExchangeChunk(comm, data.groupsets, chunk, shard);
      }
    }
    file_fluxes = {};
  }

  // Every tag must have faces somewhere, and every face sent must have arrived.
  std::vector<std::uint64_t> local_face_counts(2 * num_surfaces),
    global_face_counts(2 * num_surfaces);
  for (size_t s = 0; s < num_surfaces; ++s)
  {
    local_face_counts[2 * s] = data.surfaces[surface_tags[s]].GetNumFaces();
    local_face_counts[2 * s + 1] = faces_sent[s];
  }
  comm.all_reduce(local_face_counts.data(),
                  static_cast<int>(local_face_counts.size()),
                  global_face_counts.data(),
                  mpi::op::sum<std::uint64_t>());
  for (size_t s = 0; s < num_surfaces; ++s)
  {
    OpenSnLogicalErrorIf(global_face_counts[2 * s] == 0,
                         operation + " failed: surface " + surface_tags[s] +
                           " has no faces in any of " + std::to_string(num_partitions) + " files.");
    OpenSnLogicalErrorIf(global_face_counts[2 * s] != global_face_counts[2 * s + 1],
                         operation + " failed: " + std::to_string(global_face_counts[2 * s + 1]) +
                           " faces of " + surface_tags[s] + " were sent, but " +
                           std::to_string(global_face_counts[2 * s]) + " arrived.");
  }

  for (const auto& surface_tag : surface_tags)
  {
    auto& shard = data.surfaces[surface_tag];
    std::sort(shard.index.begin(), shard.index.end());
    for (size_t i = 1; i < shard.index.size() and error.empty(); ++i)
      if (shard.index[i].first == shard.index[i - 1].first)
        error = "surface " + surface_tag +
                " contains two faces with the same centroid key. Face centroids must differ by "
                "more than the centroid spacing " +
                std::to_string(data.centroid_spacing) + ".";
    ThrowIfAnyRankFailed(comm, error, operation);
  }

  std::uint64_t total_faces = 0;
  for (size_t s = 0; s < num_surfaces; ++s)
    total_faces += global_face_counts[2 * s];
  log.Log0Verbose1() << "Read " << total_faces << " surface faces of " << num_surfaces
                     << " surface(s) from " << num_partitions << " file(s) " << file_base
                     << "<rank>.h5 onto " << num_ranks << " rank(s).";
  return data;
}

namespace
{

/// Indices of the directions of `omegas`, three components per direction, sorted by z component.
std::vector<std::uint64_t>
SortDirectionsByZ(const std::vector<double>& omegas)
{
  std::vector<std::uint64_t> order(omegas.size() / 3);
  for (std::uint64_t i = 0; i < order.size(); ++i)
    order[i] = i;
  std::sort(order.begin(),
            order.end(),
            [&omegas](const std::uint64_t a, const std::uint64_t b)
            { return omegas[3 * a + 2] < omegas[3 * b + 2]; });
  return order;
}

/**
 * Return the direction of `omegas` nearest to `query` in Euclidean distance, skipping `exclude`,
 * and that distance. `order` is `SortDirectionsByZ(omegas)`; the search scans outward from the z
 * component of `query` and stops once the z difference alone exceeds the best distance.
 */
std::pair<std::uint64_t, double>
NearestDirection(const std::vector<double>& omegas,
                 const std::vector<std::uint64_t>& order,
                 const double* query,
                 const std::uint64_t exclude)
{
  std::pair<std::uint64_t, double> best = {exclude, std::numeric_limits<double>::infinity()};
  const auto Visit = [&](const std::uint64_t j)
  {
    if (j == exclude)
      return;
    const double* w = &omegas[3 * j];
    const double distance =
      std::sqrt((w[0] - query[0]) * (w[0] - query[0]) + (w[1] - query[1]) * (w[1] - query[1]) +
                (w[2] - query[2]) * (w[2] - query[2]));
    if (distance < best.second)
      best = {j, distance};
  };
  const auto start =
    static_cast<size_t>(std::lower_bound(order.begin(),
                                         order.end(),
                                         query[2],
                                         [&omegas](const std::uint64_t j, const double z)
                                         { return omegas[3 * j + 2] < z; }) -
                        order.begin());
  for (size_t k = start; k < order.size() and omegas[3 * order[k] + 2] - query[2] <= best.second;
       ++k)
    Visit(order[k]);
  for (size_t k = start; k > 0 and query[2] - omegas[3 * order[k - 1] + 2] <= best.second; --k)
    Visit(order[k - 1]);
  return best;
}

/**
 * Distance from each direction of `omegas` to its nearest other direction, capped at 2, the largest
 * distance between unit vectors. Throws if two directions coincide.
 */
std::vector<double>
DirectionSpacing(const std::vector<double>& omegas,
                 const std::vector<std::uint64_t>& order,
                 const std::string& description)
{
  std::vector<double> spacing(order.size());
  for (std::uint64_t j = 0; j < spacing.size(); ++j)
  {
    spacing[j] = std::min(2.0, NearestDirection(omegas, order, &omegas[3 * j], j).second);
    OpenSnInvalidArgumentIf(not(spacing[j] > 0.0),
                            "The " + description + " contains duplicate directions.");
  }
  return spacing;
}

/**
 * Map each direction of `from` to the direction of `to` with the same direction vector. Paired
 * directions must have the same weight.
 *
 * Forward direction `d` maps to its nearest adjoint direction `j` when |omega_d - omega_j| is at
 * most direction_match_tolerance times the distance from `j` to its nearest other adjoint
 * direction. Every other adjoint direction is then nearly that full spacing from `omega_d`, so the
 * match is unambiguous however closely the quadrature places its directions.
 */
std::vector<std::uint64_t>
BuildSameDirectionMap(const SurfaceQuadrature& from,
                      const SurfaceQuadrature& to,
                      const std::string& context)
{
  OpenSnInvalidArgumentIf(from.num_dirs != to.num_dirs or from.num_groups != to.num_groups or
                            from.first_group != to.first_group,
                          "The forward and adjoint data have different direction counts or group "
                          "ranges in " +
                            context + ".");
  for (const auto* omegas : {&from.omegas, &to.omegas})
    OpenSnInvalidArgumentIf(
      omegas->size() != 3 * from.num_dirs or
        not std::all_of(
          omegas->begin(), omegas->end(), [](const double value) { return std::isfinite(value); }),
      "The forward or adjoint directions of " + context + " are missing or not finite.");
  OpenSnInvalidArgumentIf(from.weights.size() != from.num_dirs or to.weights.size() != to.num_dirs,
                          "The forward or adjoint weights of " + context + " are missing.");

  const auto from_order = SortDirectionsByZ(from.omegas);
  const auto to_order = SortDirectionsByZ(to.omegas);
  DirectionSpacing(from.omegas, from_order, "forward quadrature of " + context);
  const auto to_spacing = DirectionSpacing(to.omegas, to_order, "adjoint quadrature of " + context);

  std::vector<std::uint64_t> map(from.num_dirs);
  std::vector<bool> used(to.num_dirs, false);
  for (std::uint64_t d = 0; d < from.num_dirs; ++d)
  {
    const auto [match, distance] =
      NearestDirection(to.omegas, to_order, &from.omegas[3 * d], to.num_dirs);
    if (not(distance <= direction_match_tolerance * to_spacing[match]))
    {
      std::ostringstream message;
      message << std::setprecision(3) << "Forward direction " << d << " of " << context
              << " has no matching adjoint direction. The nearest, adjoint direction " << match
              << ", is " << distance << " away, more than " << direction_match_tolerance
              << " times its spacing " << to_spacing[match] << " to the next adjoint direction.";
      OpenSnInvalidArgument(message.str());
    }
    OpenSnInvalidArgumentIf(used[match],
                            "The direction mapping of " + context + " is not one-to-one.");
    OpenSnInvalidArgumentIf(std::abs(from.weights[d] - to.weights[match]) > same_problem_tolerance,
                            "Forward direction " + std::to_string(d) + " of " + context +
                              " and its matching adjoint direction have different weights.");
    used[match] = true;
    map[d] = match;
  }
  return map;
}

/**
 * Pair the nodes of face `anchor_face` of `anchor` with those of face `partner_face` of `partner`
 * by location, within 1e-6 times the face size. `node_map[i]` is the partner node of anchor node
 * `i`.
 */
bool
MatchFaceNodes(const SurfaceShard& anchor,
               const std::uint64_t anchor_face,
               const SurfaceShard& partner,
               const std::uint64_t partner_face,
               std::vector<std::uint64_t>& node_map)
{
  const auto num_face_nodes = anchor.GetNumFaceNodes(anchor_face);
  if (partner.GetNumFaceNodes(partner_face) != num_face_nodes)
    return false;
  node_map.assign(num_face_nodes, 0);
  // A single-node face is identified by its key alone.
  if (num_face_nodes == 1)
    return true;

  const double* anchor_coords = &anchor.coords[3 * anchor.node_offset[anchor_face]];
  const double* partner_coords = &partner.coords[3 * partner.node_offset[partner_face]];
  const auto DistanceSquared = [](const double* x, const double* y)
  {
    return (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
           (x[2] - y[2]) * (x[2] - y[2]);
  };
  double face_size_sq = 0.0;
  for (std::uint64_t i = 0; i < num_face_nodes; ++i)
    for (std::uint64_t j = i + 1; j < num_face_nodes; ++j)
      face_size_sq =
        std::max(face_size_sq, DistanceSquared(&anchor_coords[3 * i], &anchor_coords[3 * j]));
  const double tolerance_sq = 1.0e-12 * face_size_sq;

  std::vector<bool> used(num_face_nodes, false);
  for (std::uint64_t i = 0; i < num_face_nodes; ++i)
  {
    bool found = false;
    for (std::uint64_t j = 0; j < num_face_nodes and not found; ++j)
      if (not used[j] and
          DistanceSquared(&anchor_coords[3 * i], &partner_coords[3 * j]) <= tolerance_sq)
      {
        node_map[i] = j;
        used[j] = true;
        found = true;
      }
    if (not found)
      return false;
  }
  return true;
}

/**
 * Return the face integral of the forward times the adjoint angular flux for one direction pair,
 * summed over groups: sum_ij M_ij sum_g psi(f_i, d, g) psi^dagger(a_j, d^dagger, g). M is the
 * anchor face's mass matrix, in anchor node order. f_i = `forward_node_map[i]` (or i when it is
 * null) and a_j = `adjoint_node_map[j]`. Both fluxes are ordered by node, direction, and group.
 */
double
FaceFluxInnerProduct(const double* mass,
                     const std::uint64_t num_face_nodes,
                     const double* forward_psi,
                     const std::vector<std::uint64_t>* forward_node_map,
                     const std::uint64_t forward_dir,
                     const double* adjoint_psi,
                     const std::vector<std::uint64_t>& adjoint_node_map,
                     const std::uint64_t adjoint_dir,
                     const std::uint64_t num_dirs,
                     const std::uint64_t num_groups)
{
  double sum = 0.0;
  for (std::uint64_t i = 0; i < num_face_nodes; ++i)
  {
    const auto forward_node = forward_node_map ? (*forward_node_map)[i] : i;
    const double* forward_node_psi =
      forward_psi + (forward_node * num_dirs + forward_dir) * num_groups;
    for (std::uint64_t j = 0; j < num_face_nodes; ++j)
    {
      const double* adjoint_node_psi =
        adjoint_psi + (adjoint_node_map[j] * num_dirs + adjoint_dir) * num_groups;
      double group_sum = 0.0;
      for (std::uint64_t g = 0; g < num_groups; ++g)
        group_sum += forward_node_psi[g] * adjoint_node_psi[g];
      sum += mass[i * num_face_nodes + j] * group_sum;
    }
  }
  return sum;
}

/// Pointer to the angular flux of face `face` of `shard` in groupset `gs`.
const double*
FacePsi(const SurfaceShard& shard,
        const std::uint64_t face,
        const size_t gs,
        const SurfaceQuadrature& quadrature)
{
  return &shard.groupsets[gs]
            .psi[shard.node_offset[face] * quadrature.num_dirs * quadrature.num_groups];
}

bool
EndsWith(const std::string& value, const std::string& suffix)
{
  return value.size() > suffix.size() and value.ends_with(suffix);
}

/// The tag of the other side of interior-surface tag `tag`, which ends in `_u` or `_d`.
std::string
OtherSideTag(const std::string& tag)
{
  return tag.substr(0, tag.size() - 2) + (EndsWith(tag, "_u") ? "_d" : "_u");
}

} // namespace

std::vector<SurfaceResponse>
ComputeSurfaceResponses(const DistributedSurfaceData& forward,
                        const DistributedSurfaceData& adjoint,
                        const std::vector<std::string>& surfaces,
                        const mpi::Communicator& comm)
{
  OpenSnInvalidArgumentIf(forward.surfaces.empty() or adjoint.surfaces.empty(),
                          "Both data sets must contain surface angular fluxes.");
  OpenSnInvalidArgumentIf(forward.adjoint,
                          "The forward data were written from a problem in adjoint mode.");
  OpenSnInvalidArgumentIf(not adjoint.adjoint,
                          "The adjoint data were written from a problem in forward mode.");
  OpenSnInvalidArgumentIf(forward.centroid_origin != adjoint.centroid_origin or
                            forward.centroid_spacing != adjoint.centroid_spacing,
                          "The forward and adjoint data have different face-centroid origins or "
                          "spacings. Both must be written for the same mesh.");
  OpenSnInvalidArgumentIf(forward.groupsets.size() != adjoint.groupsets.size(),
                          "The forward and adjoint data have different numbers of groupsets.");

  const auto num_groupsets = forward.groupsets.size();
  std::vector<std::vector<std::uint64_t>> direction_maps;
  direction_maps.reserve(num_groupsets);
  for (size_t gs = 0; gs < num_groupsets; ++gs)
    direction_maps.push_back(BuildSameDirectionMap(
      forward.groupsets[gs], adjoint.groupsets[gs], "groupset " + std::to_string(gs)));

  // Map the available surfaces by name from the forward data. The writer records whether each tag
  // is one side, <name>_u or <name>_d, of an interior surface; both sides form the interior surface
  // <name>, which DistributeSurfaceAngularFluxes reads together. Any other tag is a boundary,
  // whatever its name.
  // ex; available_surfs = { "inter_x" -> interior, "xmax" -> boundary }
  std::map<std::string, SurfaceKind> available_surfs;
  for (const auto& surface : forward.surfaces)
  {
    const auto& surface_tag = surface.first;
    std::string name = surface_tag;
    SurfaceKind kind = SurfaceKind::boundary;
    if (surface.second.interior)
    {
      const auto other = forward.surfaces.find(OtherSideTag(surface_tag));
      OpenSnInvalidArgumentIf(other == forward.surfaces.end() or not other->second.interior,
                              "The forward data hold side " + surface_tag +
                                " of an interior surface without its other side.");
      name = surface_tag.substr(0, surface_tag.size() - 2);
      kind = SurfaceKind::interior;
    }
    const auto [it, inserted] = available_surfs.emplace(name, kind);
    OpenSnInvalidArgumentIf(not inserted and it->second != kind,
                            "Surface name " + name +
                              " refers to both a boundary and an interior surface.");
  }

  // A surface to evaluate: the name it was requested by, its kind, the name of its stored tags, and
  // whether an interior surface is measured along its `_d` normal instead of its `_u` normal.
  // ex; "inter_x_d" -> { "inter_x_d", interior, "inter_x", true }
  struct SelectedSurface
  {
    std::string name;
    SurfaceKind kind = SurfaceKind::boundary;
    std::string stored_name;
    bool reversed = false;
  };
  std::vector<SelectedSurface> selected_surfs;
  if (surfaces.empty())
    // If no surfaces are specified all are used, interior surfaces along their `_u` normals.
    for (const auto& [name, kind] : available_surfs)
      selected_surfs.push_back({name, kind, name, false});
  else
  {
    // Guard for when a surface is requested multiple times, under any of its names.
    std::set<std::string> requested_surfs;
    for (const auto& surface : surfaces)
    {
      SelectedSurface selected{surface, SurfaceKind::boundary, surface, false};
      auto it = available_surfs.find(surface);
      // A side tag that is not itself an available surface names the interior surface of that
      // side, measured along that side's normal.
      if (it == available_surfs.end() and (EndsWith(surface, "_u") or EndsWith(surface, "_d")))
      {
        selected.stored_name = surface.substr(0, surface.size() - 2);
        selected.reversed = EndsWith(surface, "_d");
        it = available_surfs.find(selected.stored_name);
        if (it != available_surfs.end() and it->second != SurfaceKind::interior)
          it = available_surfs.end();
      }
      OpenSnInvalidArgumentIf(it == available_surfs.end(),
                              "Surface " + surface + " is not in the forward data.");
      selected.kind = it->second;
      OpenSnInvalidArgumentIf(not requested_surfs.insert(selected.stored_name).second,
                              "Surface " + selected.stored_name + " was requested more than once.");
      selected_surfs.push_back(std::move(selected));
    }
  }

  // The adjoint tags each surface pairs with: both sides of an interior surface, and the boundary
  // itself.
  const auto AdjointTags = [](const std::string& name, const SurfaceKind kind)
  {
    if (kind == SurfaceKind::interior)
      return std::vector<std::string>{name + "_u", name + "_d"};
    return std::vector<std::string>{name};
  };
  for (const auto& selected : selected_surfs)
    for (const auto& surface_tag : AdjointTags(selected.stored_name, selected.kind))
    {
      const auto partner = adjoint.surfaces.find(surface_tag);
      OpenSnInvalidArgumentIf(partner == adjoint.surfaces.end(),
                              std::string("Surface ")
                                .append(selected.name)
                                .append(" requires tag ")
                                .append(surface_tag)
                                .append(" in the adjoint data."));
      OpenSnInvalidArgumentIf(
        partner->second.interior != (selected.kind != SurfaceKind::boundary),
        std::string("Tag ")
          .append(surface_tag)
          .append(" of the adjoint data is ")
          .append(partner->second.interior ? "an interior surface side" : "a boundary")
          .append(", but surface ")
          .append(selected.name)
          .append(" of the forward data is not."));
    }

  // Per surface: outgoing and incoming contributions, and counts of anchor faces without a partner,
  // partner faces without an anchor, and matched faces with different node locations.
  const auto num_surfaces = selected_surfs.size();
  std::vector<double> local_responses(2 * num_surfaces, 0.0);
  std::vector<std::uint64_t> local_counts(3 * num_surfaces, 0);
  std::vector<std::uint64_t> partner_node_map, fwd_d_node_map, adj_u_node_map, adj_d_node_map;

  for (size_t surface_index = 0; surface_index < num_surfaces; ++surface_index)
  {
    const auto& name = selected_surfs[surface_index].stored_name;
    const auto kind = selected_surfs[surface_index].kind;
    double& outgoing = local_responses[2 * surface_index];
    double& incoming = local_responses[2 * surface_index + 1];
    std::uint64_t& missing = local_counts[3 * surface_index];
    std::uint64_t& extra = local_counts[3 * surface_index + 1];
    std::uint64_t& mismatched = local_counts[3 * surface_index + 2];

    // A boundary: the forward faces are the anchor and are summed against the adjoint faces of the
    // same tag. Both hold upwind values, so the outflow directions pair the forward cell value with
    // the adjoint boundary value and the inflow directions the forward boundary value with the
    // adjoint cell value.
    if (kind == SurfaceKind::boundary)
    {
      const auto& anchor = forward.surfaces.at(name);
      const auto& partner = adjoint.surfaces.at(AdjointTags(name, kind).front());
      std::uint64_t num_matched = 0;
      for (const auto& [key, anchor_face] : anchor.index)
      {
        const auto partner_face = partner.Find(key);
        if (not partner_face.has_value())
        {
          ++missing;
          continue;
        }
        ++num_matched;
        if (not MatchFaceNodes(anchor, anchor_face, partner, *partner_face, partner_node_map))
        {
          ++mismatched;
          continue;
        }

        const auto num_face_nodes = anchor.GetNumFaceNodes(anchor_face);
        const double* mass = &anchor.mass[anchor.mass_offset[anchor_face]];
        for (size_t gs = 0; gs < num_groupsets; ++gs)
        {
          const auto& quadrature = forward.groupsets[gs];
          const auto num_dirs = quadrature.num_dirs;
          const double* mu = &anchor.groupsets[gs].mu[anchor_face * num_dirs];
          const double* forward_psi = FacePsi(anchor, anchor_face, gs, quadrature);
          const double* adjoint_psi = FacePsi(partner, *partner_face, gs, adjoint.groupsets[gs]);
          for (std::uint64_t dir = 0; dir < num_dirs; ++dir)
          {
            // mu_d w_d: the cosine with the anchor normal times the direction weight.
            const auto direction_factor = mu[dir] * quadrature.weights[dir];
            if (mu[dir] == 0.0)
              continue;
            (mu[dir] > 0.0 ? outgoing : incoming) +=
              direction_factor * FaceFluxInnerProduct(mass,
                                                      num_face_nodes,
                                                      forward_psi,
                                                      nullptr,
                                                      dir,
                                                      adjoint_psi,
                                                      partner_node_map,
                                                      direction_maps[gs][dir],
                                                      num_dirs,
                                                      quadrature.num_groups);
          }
        }
      }
      extra += partner.GetNumFaces() - num_matched;
      continue;
    }

    // The anchor is the forward `_u` face; its normal is the positive axis.
    const auto& fwd_u = forward.surfaces.at(name + "_u");
    const auto& fwd_d = forward.surfaces.at(name + "_d");
    const auto& adj_u = adjoint.surfaces.at(name + "_u");
    const auto& adj_d = adjoint.surfaces.at(name + "_d");
    // Faces of fwd_d, adj_u, and adj_d matched to a fwd_u face.
    std::array<std::uint64_t, 3> num_matched = {0, 0, 0};
    for (const auto& [key, fwd_u_face] : fwd_u.index)
    {
      const auto fwd_d_face = fwd_d.Find(key);
      const auto adj_u_face = adj_u.Find(key);
      const auto adj_d_face = adj_d.Find(key);
      num_matched[0] += fwd_d_face.has_value() ? 1 : 0;
      num_matched[1] += adj_u_face.has_value() ? 1 : 0;
      num_matched[2] += adj_d_face.has_value() ? 1 : 0;
      if (not(fwd_d_face.has_value() and adj_u_face.has_value() and adj_d_face.has_value()))
      {
        ++missing;
        continue;
      }
      if (not MatchFaceNodes(fwd_u, fwd_u_face, fwd_d, *fwd_d_face, fwd_d_node_map) or
          not MatchFaceNodes(fwd_u, fwd_u_face, adj_u, *adj_u_face, adj_u_node_map) or
          not MatchFaceNodes(fwd_u, fwd_u_face, adj_d, *adj_d_face, adj_d_node_map))
      {
        ++mismatched;
        continue;
      }

      const auto num_face_nodes = fwd_u.GetNumFaceNodes(fwd_u_face);
      const double* mass = &fwd_u.mass[fwd_u.mass_offset[fwd_u_face]];
      for (size_t gs = 0; gs < num_groupsets; ++gs)
      {
        const auto& quadrature = forward.groupsets[gs];
        const auto& adj_quadrature = adjoint.groupsets[gs];
        const auto num_dirs = quadrature.num_dirs;
        const auto num_groups = quadrature.num_groups;
        const double* mu = &fwd_u.groupsets[gs].mu[fwd_u_face * num_dirs];
        const double* fwd_u_psi = FacePsi(fwd_u, fwd_u_face, gs, quadrature);
        const double* fwd_d_psi = FacePsi(fwd_d, *fwd_d_face, gs, quadrature);
        const double* adj_u_psi = FacePsi(adj_u, *adj_u_face, gs, adj_quadrature);
        const double* adj_d_psi = FacePsi(adj_d, *adj_d_face, gs, adj_quadrature);
        for (std::uint64_t dir = 0; dir < num_dirs; ++dir)
        {
          // mu_d w_d: the cosine with the anchor normal times the direction weight.
          const auto direction_factor = mu[dir] * quadrature.weights[dir];
          const auto adjoint_dir = direction_maps[gs][dir];
          // Forward leaves the `_u` cell and the adjoint leaves the `_d` cell for mu > 0, and the
          // reverse for mu < 0.
          if (mu[dir] > 0.0)
            outgoing += direction_factor * FaceFluxInnerProduct(mass,
                                                                num_face_nodes,
                                                                fwd_u_psi,
                                                                nullptr,
                                                                dir,
                                                                adj_d_psi,
                                                                adj_d_node_map,
                                                                adjoint_dir,
                                                                num_dirs,
                                                                num_groups);
          else if (mu[dir] < 0.0)
            incoming += direction_factor * FaceFluxInnerProduct(mass,
                                                                num_face_nodes,
                                                                fwd_d_psi,
                                                                &fwd_d_node_map,
                                                                dir,
                                                                adj_u_psi,
                                                                adj_u_node_map,
                                                                adjoint_dir,
                                                                num_dirs,
                                                                num_groups);
        }
      }
    }
    extra += (fwd_d.GetNumFaces() - num_matched[0]) + (adj_u.GetNumFaces() - num_matched[1]) +
             (adj_d.GetNumFaces() - num_matched[2]);
  }

  std::vector<double> global_responses(local_responses.size(), 0.0);
  std::vector<std::uint64_t> global_counts(local_counts.size(), 0);
  comm.all_reduce(local_responses.data(),
                  static_cast<int>(local_responses.size()),
                  global_responses.data(),
                  mpi::op::sum<double>());
  comm.all_reduce(local_counts.data(),
                  static_cast<int>(local_counts.size()),
                  global_counts.data(),
                  mpi::op::sum<std::uint64_t>());

  // The counts are global, so every rank reaches the same decision.
  std::string message;
  for (size_t surface_index = 0; surface_index < num_surfaces; ++surface_index)
  {
    const auto* counts = &global_counts[3 * surface_index];
    if (counts[0] == 0 and counts[1] == 0 and counts[2] == 0)
      continue;
    message += " Surface " + selected_surfs[surface_index].name + ": " + std::to_string(counts[0]) +
               " forward faces without a matching face, " + std::to_string(counts[1]) +
               " faces without a matching forward face, and " + std::to_string(counts[2]) +
               " matched faces with different node locations.";
  }
  OpenSnLogicalErrorIf(not message.empty(),
                       "The forward and adjoint surface data do not cover the same faces." +
                         message);

  // Along the `_d` normal every cosine changes sign, so the directions leaving the `_d` cell become
  // the outgoing ones and both contributions change sign.
  std::vector<SurfaceResponse> surface_responses;
  surface_responses.reserve(num_surfaces);
  for (size_t surface_index = 0; surface_index < num_surfaces; ++surface_index)
  {
    const auto& selected = selected_surfs[surface_index];
    const auto outgoing = global_responses[2 * surface_index];
    const auto incoming = global_responses[2 * surface_index + 1];
    if (selected.reversed)
      surface_responses.push_back({selected.name, selected.kind, -incoming, -outgoing});
    else
      surface_responses.push_back({selected.name, selected.kind, outgoing, incoming});
  }
  return surface_responses;
}

} // namespace opensn
